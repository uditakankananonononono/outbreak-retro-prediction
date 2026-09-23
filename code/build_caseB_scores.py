#!/usr/bin/env python3
"""Case B predictor-side: pre-2010 H3N2 HA MSA + per-site signals in H3 mature numbering.
Outputs data/predictor/caseB_site_signals.{json,csv}. Predictor-side only."""
import json, math, re, subprocess, time, urllib.request
from Bio import SeqIO, Align, Seq
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, HSExposureCB
import numpy as np

OUT = "data/predictor"
os_join = __import__("os").path.join

def get(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "open-data-research/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                time.sleep(0.4); return r.read()
        except Exception as ex:
            print(" retry", ex, flush=True); time.sleep(2)
    raise RuntimeError(url)

# --- A/Aichi/2/1968 HA (numbering reference, 1968 = predictor-side) ---
q = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&retmode=json&retmax=5&term="
     "Influenza%20A%20virus%5BOrganism%5D%20AND%20%22A%2FAichi%2F2%2F1968%22%5BAll%20Fields%5D%20AND%20hemagglutinin%5BAll%20Fields%5D%20AND%201600%3A1800%5BSLEN%5D")
ids = json.loads(get(q))["esearchresult"]["idlist"]
fa = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={ids[0]}&rettype=fasta_cds_aa&retmode=text").decode()
open(f"{OUT}/aichi68_ha.fasta", "w").write(fa)
ref_rec = next(SeqIO.parse(f"{OUT}/aichi68_ha.fasta", "fasta"))
aichi_full = str(ref_rec.seq)
# mature numbering: strip 16-aa signal peptide; HA1 then starts with QDLP...
gi = aichi_full.find("GNDNSTAT")
assert gi > 3, "signal peptide cleavage not found"
idx = gi - 4
aichi_mature = aichi_full[idx:]  # position 1 = Q (H3 mature numbering)
print("Aichi68 HA mature start at offset", idx, "len", len(aichi_mature))

# --- ORF extraction + year-stratified subsample ---
meta = json.load(open(f"{OUT}/h3n2_metadata.json"))
acc2meta = {m["accession"]: m for m in meta.values()}
by_year = {}
kept_orf = {}
for rec in SeqIO.parse(f"{OUT}/h3n2_ha.fasta", "fasta"):
    acc = rec.id.split("|")[0]
    nt = str(rec.seq).upper().replace("U", "T")
    best = None
    for frame in range(3):
        aa = str(Seq.Seq(nt[frame: len(nt) - (len(nt) - frame) % 3]).translate())
        # HA CDS ~566 aa; find window without stops
        cand = aa[:aa.find("*")] if "*" in aa else aa
        if best is None or len(cand) > len(best): best = cand
    if best is None or len(best) < 540: continue
    kept_orf[acc] = best[:570]
    t = acc2meta.get(acc, {}).get("title", "")
    m = re.search(r"/((?:19|20)\d{2})", t)
    yr = m.group(1) if m else acc2meta.get(acc, {}).get("createdate", "9999")[:4]
    by_year.setdefault(yr, []).append(acc)
print("ORF kept:", len(kept_orf), "years:", len(by_year))
sub = []
for yr in sorted(by_year):
    sub.extend(by_year[yr][:40])
print("subsample:", len(sub))

# --- MSA ---
with open("work/caseB_ha.fasta", "w") as fh:
    fh.write(">AICHI68_ref\n" + aichi_mature + "\n")
    for acc in sub:
        fh.write(f">{acc}\n{kept_orf[acc]}\n")
subprocess.run(["/home/sandbox/bin/muscle", "-super5", "work/caseB_ha.fasta",
                "-output", "work/caseB_ha_msa.fasta", "-threads", "2"], check=True)
msa = list(SeqIO.parse("work/caseB_ha_msa.fasta", "fasta"))
ref = next(r for r in msa if r.id == "AICHI68_ref")
others = [str(r.seq) for r in msa if r.id != "AICHI68_ref"]
print("MSA:", len(others), "x", len(ref.seq))
pos2col, p = {}, 0
for j, c in enumerate(str(ref.seq)):
    if c != "-":
        p += 1; pos2col[p] = j
HA1_END = 329

s1, cov = {}, {}
from collections import Counter
for pos in range(1, HA1_END + 1):
    j = pos2col[pos]
    chars = [s[j] for s in others if s[j] != "-"]
    cov[pos] = len(chars) / len(others)
    if len(chars) >= 10:
        cnt = Counter(chars); n = len(chars)
        s1[pos] = -sum((v/n) * math.log(v/n) for v in cnt.values())

# --- pre-2010 H3 structure ---
cands = ["1HGD", "2HMG", "1HGF", "3HMG", "4HMG", "1EO8", "2VIR", "1KEN"]
pdb = None
for c in cands:
    info = json.loads(get(f"https://data.rcsb.org/rest/v1/core/entry/{c}"))
    rel = info.get("rcsb_accession_info", {}).get("initial_release_date", "9999")
    title = info.get("struct", {}).get("title", "").upper()
    if rel <= "2009-12-31" and ("HEMAGGLUTININ" in title or "HAEMAGGLUTININ" in title):
        pdb = c
        print("structure pick:", c, rel, title[:70])
        open(f"{OUT}/{c}.cif", "wb").write(get(f"https://files.rcsb.org/download/{c}.cif"))
        break
assert pdb, "no pre-2010 H3 HA structure found"
parser = MMCIFParser(QUIET=True)
st = parser.get_structure(pdb, f"{OUT}/{pdb}.cif")
model = st[0]
hse = HSExposureCB(model)
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -8; aligner.extend_gap_score = -0.5; aligner.mode = "local"
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
best = None
for ch in model:
    resids = [r.id[1] for r in ch if r.id[0] == " " and is_aa(r)]
    cseq = "".join(seq1(r.resname) for r in ch if r.id[0] == " " and is_aa(r))
    if len(cseq) < 200: continue
    aln = aligner.align(aichi_mature[:HA1_END], cseq)[0]
    m = {}
    rb, tb = aln.aligned[0], aln.aligned[1]
    for (rs, re_), (ts, te) in zip(rb, tb):
        for k in range(re_ - rs): m[int(resids[ts + k])] = int(rs + k + 1)
    if best is None or len(m) > len(best[0]): best = (m, ch.id)
mapS, chS = best
print(pdb, "chain", chS, "mapped HA1 residues:", len(mapS))
hse_dict = {(r.get_parent().id, r.id[1]): r.xtra.get("EXP_HSE_B_U") for r in model.get_residues()}
s2 = {int(pos): float(hse_dict.get((chS, resid))) for resid, pos in mapS.items()
      if pos <= HA1_END and hse_dict.get((chS, resid)) is not None}

# s3: distance to RBS footprint (pre-2010 literature: 98,153,183,190,194,195,226,228)
foot = {98, 153, 183, 190, 194, 195, 226, 228}
pos2resid = {p: r for r, p in mapS.items()}
foot_coords = []
for fp in foot:
    if fp in pos2resid:
        try: foot_coords.append(model[chS][pos2resid[fp]]["CA"].coord)
        except KeyError: pass
s3 = {}
for pos, resid in mapS.items():
    if pos > HA1_END: continue
    try: ca = model[chS][resid]["CA"].coord
    except KeyError: continue
    s3[int(pos)] = -min(float(np.linalg.norm(ca - fc)) for fc in foot_coords)

# s4: epitope membership, Wilson & Cox 1990 ranges (pre-2010), H3 numbering
EPIT = {"A": [(121,147)], "B": [(155,160),(186,198)], "C": [(44,54),(273,278)],
        "D": [(96,98),(201,227),(243,244)], "E": [(26,27),(57,63),(75,83),(91,92)]}
s4 = {}
for pos in range(1, HA1_END + 1):
    s4[pos] = 1.0 if any(lo <= pos <= hi for rs in EPIT.values() for lo, hi in rs) else 0.0

json.dump({"s1_entropy": s1, "s2_hse_up": s2, "s3_neg_rbs_dist": s3, "s4_epitope": s4,
           "coverage": cov, "structure": pdb, "n_msa": len(others)},
          open(f"{OUT}/caseB_site_signals.json", "w"))
import csv
with open(f"{OUT}/caseB_site_signals.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["h3_pos", "coverage", "s1_entropy", "s2_hse_up", "s3_neg_rbs_dist", "s4_epitope"])
    for pos in range(1, HA1_END + 1):
        w.writerow([pos, round(cov.get(pos, 0), 3), s1.get(pos), s2.get(pos), s3.get(pos), s4.get(pos)])
print("caseB signals written")
