#!/usr/bin/env python3
"""T2 variant-interpret: annotate a variant's mutations with pre-registered risk context.
Usage: python3 tools/variant_interpret.py --case sarscov2 --mutations S:E484K,S:N501Y,S:K417N
       python3 tools/variant_interpret.py --case h3n2 --mutations HA1:N145D,HA1:K160T
Output: TSV table. Risk tier from the validated PROMIS score percentile (pre-outbreak predictors)."""
import argparse, json, re, os

ap = argparse.ArgumentParser()
ap.add_argument("--case", choices=["sarscov2", "h3n2"], required=True)
ap.add_argument("--mutations", required=True)
a = ap.parse_args()
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
comp = json.load(open(f"{root}/data/predictor/composite_scores.json"))

if a.case == "sarscov2":
    compA = {int(k): v for k, v in comp["composite_A_sars1"].items()}
    m12 = {int(k): int(v) for k, v in json.load(open(f"{root}/results/sars1_to_sars2_spike_map.json")).items()}
    m21 = {v: k for k, v in m12.items()}
    vals = sorted(compA.values())
    def lookup(pos):
        s1p = m21.get(pos)
        if s1p is None or s1p not in compA: return None, None
        z = compA[s1p]
        pct = 100 * sum(1 for v in vals if v <= z) / len(vals)
        return s1p, (z, pct)
else:
    compB = {int(k): v for k, v in comp["composite_B_h3"].items()}
    vals = sorted(compB.values())
    def lookup(pos):
        if pos not in compB: return None, None
        z = compB[pos]
        pct = 100 * sum(1 for v in vals if v <= z) / len(vals)
        return pos, (z, pct)

def tier(pct):
    return "HIGH" if pct >= 90 else ("WATCH" if pct >= 75 else "background")

print("mutation\tposition\tref_position\tcomposite_z\tpercentile\ttier")
for mut in a.mutations.split(","):
    m = re.match(r"(?:S|HA1):([A-Z])(\d+)([A-Z])", mut.strip())
    if not m: print(f"{mut}\t?\t?\t?\t?\tunparsed"); continue
    pos = int(m.group(2))
    rp, res = lookup(pos)
    if res is None:
        print(f"{mut}\t{pos}\t{rp or '-'}\t-\t-\tno-data")
    else:
        z, pct = res
        print(f"{mut}\t{pos}\t{rp}\t{z:.2f}\t{pct:.0f}%\t{tier(pct)}")
