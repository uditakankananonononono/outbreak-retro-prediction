#!/usr/bin/env python3
"""Composite scores + PREREGISTRATION.md (hashed). Predictor-side only. Run BEFORE outcome data."""
import json, hashlib, math, datetime

def zmap(d, invert=False):
    vals = [(p, (-v if invert else v)) for p, v in d.items() if v is not None]
    xs = [v for _, v in vals]
    mu = sum(xs) / len(xs); sd = (sum((x - mu) ** 2 for x in xs) / len(xs)) ** 0.5 or 1.0
    return {int(p): (v - mu) / sd for p, v in vals}

A = json.load(open("data/predictor/caseA_site_signals.json"))
zA = {"s1": zmap(A["s1_entropy"]), "s2": zmap(A["s2_hse_up"], invert=True),
      "s3": zmap(A["s3_neg_ace2_dist"]), "s4": zmap(A["s4_hostshift"])}
covA = {int(k): v for k, v in A["coverage"].items()}
compA = {}
for pos in range(1, A["spike_len"] + 1):
    if covA.get(pos, 0) < 0.5: continue
    zs = [zA[s][pos] for s in ("s1", "s2", "s3", "s4") if pos in zA[s]]
    if len(zs) >= 3:
        compA[pos] = sum(zs) / len(zs)
rbd = sorted(((p, s) for p, s in compA.items() if 318 <= p <= 510), key=lambda t: -t[1])
ntd = sorted(((p, s) for p, s in compA.items() if 14 <= p <= 305), key=lambda t: -t[1])
top_rbd, top_ntd = rbd[:30], ntd[:20]

B = json.load(open("data/predictor/caseB_site_signals.json"))
zB = {"s1": zmap(B["s1_entropy"]), "s2": zmap(B["s2_hse_up"], invert=True),
      "s3": zmap(B["s3_neg_rbs_dist"]), "s4": zmap(B["s4_epitope"])}
covB = {int(k): v for k, v in B["coverage"].items()}
compB = {}
for pos in range(1, 330):
    if covB.get(pos, 0) < 0.5: continue
    zs = [zB[s][pos] for s in ("s1", "s2", "s3", "s4") if pos in zB[s]]
    if len(zs) >= 3:
        compB[pos] = sum(zs) / len(zs)
top_ha1 = sorted(compB.items(), key=lambda t: -t[1])[:30]

json.dump({"composite_A_sars1": compA, "composite_B_h3": compB},
          open("data/predictor/composite_scores.json", "w"))

lines = []
lines.append("# PREREGISTRATION - Project 13 retrospective outbreak prediction")
lines.append(f"Generated: {datetime.datetime.utcnow().isoformat()}Z. Committed BEFORE any outcome data download or scoring.")
lines.append("")
lines.append("## Method (as locked in GATES_LOCKED.md)")
lines.append("Composite site-risk score PROMIS = unweighted mean of z-scored signals s1 variability,")
lines.append("s2 exposure (HSE inverted), s3 receptor proximity (negative min-CA distance), s4 context")
lines.append("(Case A bat-vs-human consensus shift; Case B pre-2010 epitope membership). Site eligible if")
lines.append("MSA coverage >= 0.5 and >= 3 of 4 signals present. No weights fitted anywhere.")
lines.append("")
lines.append("## Scoring rules (locked)")
lines.append("- G3: Case A enrichment of VOC-defining spike substitutions (WHO VOCs Alpha, Beta, Gamma,")
lines.append("  Delta, Omicron BA.1+BA.2, from cov-lineages constellations) among predicted sites vs")
lines.append("  10,000-draw label-permutation null over all scored RBD+NTD sites; one-sided p<0.05.")
lines.append("- G4: Case B enrichment of HA1 substitutions reaching >50% global frequency (first")
lines.append("  appearance 2010-2019, Nextstrain 12y tree + tip frequencies) among predicted sites; same null.")
lines.append("- G5: site-level ROC-AUC of composite (Case A, positives = VOC sites in RBD+NTD) > AUC of")
lines.append("  every baseline B1 variability-only, B2 exposure-only, B3 Hie-style 3-mer grammaticality +")
lines.append("  BLOSUM semantic change, B4 EVEscape-style independent-site frequency x exposure;")
lines.append("  10,000-draw bootstrap 95% CI on AUC difference vs best baseline excludes 0.")
lines.append("- Coordinate mapping SARS1->SARS2 spike via pairwise alignment of the two reference")
lines.append("  sequences (NC_045512 used ONLY as scoring coordinate frame), applied at scoring time.")
lines.append("")
lines.append("## Case A predictions (SARS-CoV-1 spike numbering; mapped to SARS-CoV-2 only at scoring)")
lines.append("Top-30 RBD (318-510): " + ", ".join(str(p) for p, _ in top_rbd))
lines.append("Top-20 NTD (14-305): " + ", ".join(str(p) for p, _ in top_ntd))
lines.append("")
lines.append("## Case B predictions (H3 mature HA numbering, HA1 1-329)")
lines.append("Top-30 HA1: " + ", ".join(str(p) for p, _ in top_ha1))
lines.append("")
lines.append("## Predicted-site scores")
lines.append("RBD: " + "; ".join(f"{p}:{s:.2f}" for p, s in top_rbd))
lines.append("NTD: " + "; ".join(f"{p}:{s:.2f}" for p, s in top_ntd))
lines.append("HA1: " + "; ".join(f"{p}:{s:.2f}" for p, s in top_ha1))
body = "\n".join(lines) + "\n"
h = hashlib.sha256(body.encode()).hexdigest()
with open("PREREGISTRATION.md", "w") as fh:
    fh.write(body + f"\nSHA256(of body above) = {h}\n")
print("preregistration written; sha256", h)
print("RBD top30:", [p for p, _ in top_rbd])
print("NTD top20:", [p for p, _ in top_ntd])
print("HA1 top30:", [p for p, _ in top_ha1])
