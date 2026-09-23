#!/usr/bin/env python3
"""T1 varscan-risk: per-site pre-outbreak mutation-risk scoring.
Usage: python3 tools/varscan_risk.py --msa aligned.fasta --ref-id REFSEQID \
         [--structure file.cif --target-chain A --receptor-cif file2.cif --receptor-target-chain E --receptor-chain A] \
         --out risk.csv
Signals: s1 entropy across MSA; s2 half-sphere exposure (structure); s3 receptor proximity
(receptor complex); composite = mean of available z-scores. Pure stdlib+biopython+numpy."""
import argparse, csv, math
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument("--msa", required=True); ap.add_argument("--ref-id", required=True)
ap.add_argument("--structure"); ap.add_argument("--target-chain", default="A")
ap.add_argument("--receptor-target-chain"); ap.add_argument("--receptor-chain")
ap.add_argument("--out", required=True)
a = ap.parse_args()

from Bio import SeqIO
msa = list(SeqIO.parse(a.msa, "fasta"))
ref = next(r for r in msa if r.id == a.ref_id)
others = [str(r.seq) for r in msa if r.id != a.ref_id]
pos2col, p = {}, 0
for j, c in enumerate(str(ref.seq)):
    if c != "-": p += 1; pos2col[p] = j
L = p
s1, cov = {}, {}
for pos in range(1, L + 1):
    j = pos2col[pos]
    chars = [s[j] for s in others if s[j] != "-"]
    cov[pos] = len(chars) / max(1, len(others))
    if len(chars) >= 5:
        cnt = Counter(chars); n = len(chars)
        s1[pos] = -sum((v/n) * math.log(v/n) for v in cnt.values())
s2, s3 = {}, {}
if a.structure:
    from Bio.PDB import MMCIFParser, PDBParser, HSExposureCB
    from Bio.PDB.Polypeptide import is_aa
    from Bio.SeqUtils import seq1
    from Bio import Align
    from Bio.Align import substitution_matrices
    import numpy as np
    parser = MMCIFParser(QUIET=True) if a.structure.endswith(".cif") else PDBParser(QUIET=True)
    model = parser.get_structure("s", a.structure)[0]
    hse = HSExposureCB(model)
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -8; aligner.extend_gap_score = -0.5; aligner.mode = "local"
    ref_aa = str(ref.seq).replace("-", "")
    ch = model[a.target_chain]
    resids = [r.id[1] for r in ch if r.id[0] == " " and is_aa(r)]
    cseq = "".join(seq1(r.resname) for r in ch if r.id[0] == " " and is_aa(r))
    aln = aligner.align(ref_aa, cseq)[0]
    m = {}
    for (rs, re_), (ts, te) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(re_ - rs): m[int(resids[ts + k])] = int(rs + k + 1)
    hse_d = {(r.get_parent().id, r.id[1]): r.xtra.get("EXP_HSE_B_U") for r in model.get_residues()}
    for resid, pos in m.items():
        v = hse_d.get((a.target_chain, resid))
        if v is not None: s2[pos] = float(v)
    if a.receptor_chain:
        tch = a.receptor_target_chain or a.target_chain
        rc = [x.coord for x in model[a.receptor_chain].get_atoms() if x.name == "CA"]
        for resid, pos in m.items():
            try: ca = model[tch][resid]["CA"].coord
            except KeyError: continue
            s3[pos] = -min(float(np.linalg.norm(ca - x)) for x in rc)

def zmap(d, invert=False):
    items = [(k, (-v if invert else v)) for k, v in d.items() if v is not None]
    if not items: return {}
    xs = [v for _, v in items]; mu = sum(xs)/len(xs); sd = (sum((x-mu)**2 for x in xs)/len(xs))**0.5 or 1.0
    return {k: (v-mu)/sd for k, v in items}
z1, z2, z3 = zmap(s1), zmap(s2, invert=True), zmap(s3)
with open(a.out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["position", "coverage", "s1_variability", "s2_exposure", "s3_receptor_prox", "composite_z"])
    for pos in range(1, L + 1):
        zs = [z.get(pos) for z in (z1, z2, z3) if z.get(pos) is not None]
        comp = sum(zs)/len(zs) if zs else None
        w.writerow([pos, round(cov.get(pos, 0), 3), z1.get(pos), z2.get(pos), z3.get(pos),
                    round(comp, 4) if comp is not None else None])
print(f"varscan-risk: {L} sites scored from {len(others)} MSA sequences -> {a.out}")
