#!/usr/bin/env python3
"""Case A predictor-side: pre-2020 sarbecovirus spike MSA + per-site signals (SARS-CoV-1 numbering).
Outputs data/predictor/caseA_spike_msa.fasta and caseA_site_signals.csv. NO outcome data used."""
import json, os, re, subprocess, hashlib
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, HSExposureCB
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

OUT = "data/predictor"
os.makedirs("work", exist_ok=True)

# --- SARS-CoV-1 spike coords from the reference genbank ---
gb = SeqIO.read(f"{OUT}/sars1_ref_NC_004718.gb", "genbank")
s_cds = None
for f in gb.features:
    if f.type == "CDS" and ("S" in f.qualifiers.get("gene", []) or "spike" in str(f.qualifiers.get("product", [""])).lower()) and "surface" in str(f.qualifiers.get("product", [""])).lower():
        s_cds = f; break
if s_cds is None:
    for f in gb.features:
        if f.type == "CDS" and "spike" in str(f.qualifiers.get("product", [""])).lower():
            s_cds = f; break
assert s_cds is not None, "spike CDS not found"
ref_spike_nt = s_cds.extract(gb.seq)
ref_spike_aa = s_cds.qualifiers["translation"][0] if "translation" in s_cds.qualifiers else str(ref_spike_nt.translate())
print("SARS1 spike:", len(ref_spike_aa), "aa; CDS", s_cds.location)

# --- classify host from metadata ---
meta = json.load(open(f"{OUT}/sarbecovirus_metadata.json"))
acc2host = {}
for uid, m in meta.items():
    t = (m["title"] + " " + m["organism"]).lower()
    if "bat" in t or "rhinoloph" in t: h = "bat"
    elif "civet" in t or "raccoon" in t: h = "civet"
    elif "pangolin" in t or "manis" in t: h = "pangolin"
    elif "human" in t or "sars coronavirus" in t or "severe acute respiratory syndrome coronavirus" == m["organism"].lower(): h = "human"
    else: h = "other"
    acc2host[m["accession"]] = h

# --- extract spike from each genome by local alignment to ref spike nt ---
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("NUC.4.4")
aligner.open_gap_score = -6; aligner.extend_gap_score = -1
aligner.mode = "local"
spikes = []  # (acc, host, aa)
recs = list(SeqIO.parse(f"{OUT}/sarbecovirus_genomes.fasta", "fasta"))
print("genomes:", len(recs))
for i, rec in enumerate(recs):
    acc = rec.id.split("|")[0]
    g = str(rec.seq).upper()
    try:
        aln = aligner.align(ref_spike_nt, g)[0]
        # aligned blocks: aln.aligned[1] gives target blocks; find where ref pos 0 maps
        ref_blocks, tgt_blocks = aln.aligned[0], aln.aligned[1]
        # collect target sequence for ref-covered region
        pieces, prev = [], None
        for (rs, re_), (ts, te) in zip(ref_blocks, tgt_blocks):
            if prev is not None and ts > prev:
                pieces.append("N" * min(ts - prev, 9))  # small insertions tolerated as X
            pieces.append(g[ts:te])
            prev = te
        span = "".join(pieces)
        cov = (ref_blocks[-1][1] - ref_blocks[0][0]) / len(ref_spike_nt)
        if cov < 0.9 or len(span) < len(ref_spike_nt) * 0.9: continue
        aa = str(Seq(span[: len(span) - len(span) % 3]).translate())
        if aa.count("X") > 30 or "*" in aa[:-1]: continue
        spikes.append((acc, acc2host.get(acc, "other"), aa))
    except Exception as ex:
        pass
    if (i + 1) % 50 == 0: print(f"  extracted {len(spikes)} spikes from {i+1} genomes", flush=True)
print("spikes extracted:", len(spikes))
# dedupe identical AA
seen, uniq = {}, []
for acc, h, aa in spikes:
    if aa not in seen:
        seen[aa] = acc
        uniq.append((acc, h, aa))
print("unique spikes:", len(uniq))
json.dump({acc: h for acc, h, _ in uniq}, open(f"{OUT}/caseA_spike_hosts.json", "w"), indent=0)
with open("work/caseA_spikes.fasta", "w") as fh:
    # reference first for coordinate mapping
    fh.write(">SARS1_NC_004718_ref\n" + ref_spike_aa + "\n")
    for acc, h, aa in uniq:
        fh.write(f">{acc}|{h}\n{aa}\n")
