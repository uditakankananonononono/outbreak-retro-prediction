#!/usr/bin/env python3
"""OUTCOME-side data. Run ONLY after PREREGISTRATION.md is committed. Writes data/outcome/ + manifest."""
import json, os, time, hashlib, urllib.request

OUT = "data/outcome"
os.makedirs(OUT, exist_ok=True)

def get(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "open-data-research/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                time.sleep(0.4); return r.read()
        except Exception as ex:
            print(" retry", ex, flush=True); time.sleep(2)
    raise RuntimeError(url)

# WHO VOC constellations from cov-lineages/constellations
voc_files = {"Alpha": "cB.1.1.7.json", "Beta": "cB.1.351.json", "Gamma": "cP.1.json",
             "Delta": "cB.1.617.2.json", "Omicron_BA1": "cBA.1.json", "Omicron_BA2": "cBA.2.json"}
base = "https://raw.githubusercontent.com/cov-lineages/constellations/main/constellations/definitions/"
for name, fn in voc_files.items():
    b = get(base + fn)
    open(f"{OUT}/constellation_{name}.json", "wb").write(b)
    d = json.loads(b)
    print(name, "sites:", len(d.get("sites", [])), "lineage:", d.get("lineage_name", fn))

# SARS-CoV-2 reference (coordinate frame ONLY)
b = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=gbwithparts&retmode=text")
open(f"{OUT}/sars2_ref_NC_045512.gb", "wb").write(b)

# Nextstrain H3N2 12y tree + tip frequencies
b = get("https://data.nextstrain.org/flu_seasonal_h3n2_ha_12y.json")
open(f"{OUT}/nextstrain_h3n2_ha_12y.json", "wb").write(b)
print("nextstrain tree bytes:", len(b))
b = get("https://data.nextstrain.org/flu_seasonal_h3n2_ha_12y_tip-frequencies.json")
open(f"{OUT}/nextstrain_h3n2_tip_frequencies.json", "wb").write(b)
print("frequencies bytes:", len(b))

man = {}
for root, _, files in os.walk(OUT):
    for fn in files:
        p = os.path.join(root, fn)
        man[p] = hashlib.sha256(open(p, "rb").read()).hexdigest()
json.dump(man, open("data/outcome_manifest.json", "w"), indent=1)
print("outcome manifest:", len(man), "files")
