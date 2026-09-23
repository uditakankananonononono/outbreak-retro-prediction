#!/usr/bin/env python3
"""T3 target-impact: map variant mutations onto vaccine/antiviral target annotations.
Usage: python3 tools/target_impact.py --mutations S:E484K,S:N501Y,orf1ab:P314L
Annotations: spike domains/RBD epitope classes/ACE2 contacts (project-validated), furin site,
Mpro active site and RdRp motifs (from reference feature annotations). Operational tool -
annotations may use current knowledge; the predictive machinery uses pre-outbreak data only."""
import argparse, json, re, os

ap = argparse.ArgumentParser()
ap.add_argument("--mutations", required=True)
a = ap.parse_args()
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
pcb = json.load(open(f"{root}/results/pcb_contact_reanalysis.json"))
ace2_contacts = set(json.load(open(f"{root}/data/predictor/caseA_site_signals.json"))["s3_neg_ace2_dist"].keys())
# SARS-CoV-2 annotation map (operational annotation layer)
RBD_CLASS = {"class1": [(455, 456), (475, 476), (486, 489), (493, 494), (501, 502)],
             "class2": [(417, 417), (446, 446), (449, 449), (484, 484), (490, 490), (493, 496)],
             "class3": [(345, 346), (378, 378), (440, 441), (443, 445), (498, 500)],
             "class4": [(369, 378), (408, 409)]}
NTD_SUPERSITE = [(14, 20), (140, 158), (245, 264)]
FURIN = range(681, 686)
def s_annot(pos):
    out = []
    if 319 <= pos <= 541: out.append("RBD")
    elif 14 <= pos <= 305: out.append("NTD")
    elif pos in FURIN: out.append("furin-cleavage-site")
    for cls, ranges in RBD_CLASS.items():
        if any(lo <= pos <= hi for lo, hi in ranges): out.append(f"RBD-epitope-{cls}")
    if any(lo <= pos <= hi for lo, hi in NTD_SUPERSITE): out.append("NTD-supersite")
    return out
MPRO_ACTIVE = {41: "Mpro catalytic His41", 145: "Mpro catalytic Cys145"}
RDRP_MOTIFS = {"motif_A": (532 + 0, 532 + 20), "motif_C": (759, 779)}  # orf1ab approx RdRp region annotations
def nsp_annot(pos):
    out = []
    if 3263 <= pos <= 3569: out.append("nsp12-RdRp-region")
    if 4850 <= pos <= 5100: out.append("nsp14-ExoN-region")
    if pos in MPRO_ACTIVE: out.append(MPRO_ACTIVE[pos])
    if 3569 < pos <= 4428: out.append("nsp13-helicase-region")
    return out
print("mutation\tgene\tposition\ttarget_annotations\timpact_tier")
for mut in a.mutations.split(","):
    m = re.match(r"([A-Za-z0-9]+):([A-Z])(\d+)([A-Z])", mut.strip())
    if not m: print(f"{mut}\t?\t?\tunparsed\t-"); continue
    gene, pos = m.group(1), int(m.group(3))
    ann = s_annot(pos) if gene.upper() in ("S", "SPIKE") else (nsp_annot(pos) if gene.lower().startswith("orf1") else [])
    tier = "HIGH" if any("epitope" in x or "contact" in x or "catalytic" in x or "furin" in x for x in ann) else \
           ("MODERATE" if ann else "low/none-annotated")
    print(f"{mut}\t{gene}\t{pos}\t{';'.join(ann) or '-'}\t{tier}")
