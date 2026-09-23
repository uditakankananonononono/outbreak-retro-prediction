#!/usr/bin/env python3
"""Case A per-site signals in SARS-CoV-1 spike numbering. Predictor-side only."""
import json, math, subprocess
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, HSExposureCB
import numpy as np

OUT = "data/predictor"
# --- MSA ---
subprocess.run(["/home/sandbox/bin/muscle", "-align", "work/caseA_spikes.fasta",
                "-output", "work/caseA_spike_msa.fasta", "-threads", "2"], check=True)
msa = list(SeqIO.parse("work/caseA_spike_msa.fasta", "fasta"))
ref = next(r for r in msa if r.id == "SARS1_NC_004718_ref")
others = [r for r in msa if r.id != "SARS1_NC_004718_ref"]
print("MSA:", len(others), "seqs x", len(ref.seq), "cols")
# map: spike position (1-based) -> column
pos2col, p = {}, 0
for j, c in enumerate(str(ref.seq)):
    if c != "-":
        p += 1; pos2col[p] = j
L = p
print("spike length:", L)

hosts = json.load(open(f"{OUT}/caseA_spike_hosts.json"))
seqs = [(r.id.split("|")[0], r.description.split("|")[-1], str(r.seq)) for r in others]

def col_chars(j):
    return [s[2][j] for s in seqs]

def entropy(chars):
    chars = [c for c in chars if c != "-"]
    if len(chars) < 5: return None, len(chars)
    from collections import Counter
    cnt = Counter(chars); n = len(chars)
    return -sum((v/n) * math.log(v/n) for v in cnt.values()), n

def consensus(chars, sel):
    from collections import Counter
    c = Counter(s[2][j] for j, s in [(i, s) for i, s in enumerate(seqs)] if False)  # placeholder
    return None

# host-shift per column: bat consensus vs human consensus
def cons_at(j, hostset):
    from collections import Counter
    cs = [s[2][j] for s in seqs if s[1] in hostset and s[2][j] != "-"]
    if len(cs) < 3: return None
    return Counter(cs).most_common(1)[0][0]

s1 = {}; s4 = {}; cov = {}
for pos in range(1, L + 1):
    j = pos2col[pos]
    chars = [s[2][j] for s in seqs]
    e, n = entropy(chars)
    cov[pos] = n / len(seqs)
    if e is not None: s1[pos] = e
    cb, ch = cons_at(j, {"bat"}), cons_at(j, {"human"})
    if cb and ch: s4[pos] = 1.0 if cb != ch else 0.0

# --- structure signals ---
parser = MMCIFParser(QUIET=True)
st5 = parser.get_structure("5xlr", f"{OUT}/5XLR.cif")
model5 = st5[0]
hse = HSExposureCB(model5)
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -8; aligner.extend_gap_score = -0.5; aligner.mode = "local"
ref_aa = str(ref.seq).replace("-", "")

def map_chain_to_ref(chain, refseq):
    cseq = "".join(r.resname and __import__("Bio.PDB.Polypeptide", fromlist=["is_aa"]).is_aa(r) and __import__("Bio.SeqUtils", fromlist=["seq1"]).seq1(r.resname) or "X" for r in chain if r.id[0] == " ")
    resids = [r.id[1] for r in chain if r.id[0] == " "]
    aln = aligner.align(refseq, cseq)[0]
    rb, tb = aln.aligned[0], aln.aligned[1]
    m = {}
    for (rs, re_), (ts, te) in zip(rb, tb):
        for k in range(re_ - rs):
            m[resids[ts + k]] = rs + k + 1  # struct resid -> spike pos (1-based)
    return m, cseq

best = None
for ch in model5:
    m, cseq = map_chain_to_ref(ch, ref_aa)
    if best is None or len(m) > len(best[0]): best = (m, ch.id, cseq)
mapA, chA, _ = best
print("5XLR chain", chA, "mapped residues:", len(mapA))
s2 = {}
hse_dict = {(r.get_parent().id, r.id[1]): r.xtra.get("EXP_HSE_B_U", None) for r in model5.get_residues()}
for resid, pos in mapA.items():
    v = hse_dict.get((chA, resid))
    if v is not None: s2[pos] = float(v)  # HSE up = more buried; z-invert later

# 2AJF: find RBD chain (maps to spike 300-530) and ACE2 chain (longest)
st2 = parser.get_structure("2ajf", f"{OUT}/2AJF.cif")
model2 = st2[0]
rbd_map, ace2_coords, rbd_ch = None, None, None
rbd_ref = ref_aa[300:540]
best2 = None
for ch in model2:
    m, cseq = map_chain_to_ref(ch, ref_aa)
    hit_rbd = [p for p in m.values() if 305 <= p <= 530]
    if len(hit_rbd) > 80 and (best2 is None or len(hit_rbd) > len(best2[0])):
        best2 = (m, ch.id)
for ch in model2:
    if best2 and ch.id != best2[1]:
        ca = [a.coord for a in ch.get_atoms() if a.name == "CA"]
        if ace2_coords is None or len(ca) > len(ace2_coords): ace2_coords, _ = (ca, ch.id)
mapE, chE = best2
print("2AJF RBD chain", chE, "mapped:", len(mapE))
s3 = {}
for resid, pos in mapE.items():
    try:
        ca = model2[chE][resid]["CA"].coord
    except KeyError: continue
    d = min(float(np.linalg.norm(ca - a2)) for a2 in ace2_coords)
    s3[pos] = -d  # closer = higher risk signal

json.dump({"s1_entropy": s1, "s2_hse_up": s2, "s3_neg_ace2_dist": s3, "s4_hostshift": s4,
           "coverage": cov, "spike_len": L}, open(f"{OUT}/caseA_site_signals.json", "w"))
import csv
with open(f"{OUT}/caseA_site_signals.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["sars1_spike_pos", "domain", "coverage", "s1_entropy", "s2_hse_up", "s3_neg_ace2_dist", "s4_hostshift"])
    for pos in range(1, L + 1):
        dom = "NTD" if 14 <= pos <= 305 else ("RBD" if 318 <= pos <= 510 else ("S1other" if pos <= 667 else "S2"))
        w.writerow([pos, dom, round(cov.get(pos, 0), 3), s1.get(pos), s2.get(pos), s3.get(pos), s4.get(pos)])
print("signals written for", L, "sites")
