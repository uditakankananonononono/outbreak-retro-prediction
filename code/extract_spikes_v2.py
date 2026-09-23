#!/usr/bin/env python3
"""Extract spike AA directly from GenBank CDS annotations (replaces pairwise extraction)."""
import json, time, urllib.parse, urllib.request
from Bio import SeqIO
import io

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = "data/predictor"
meta = json.load(open(f"{OUT}/sarbecovirus_metadata.json"))
acc2host = {}
for uid, m in meta.items():
    t = (m["title"] + " " + m["organism"]).lower()
    if "bat" in t or "rhinoloph" in t: h = "bat"
    elif "civet" in t or "raccoon" in t: h = "civet"
    elif "pangolin" in t or "manis" in t: h = "pangolin"
    elif "human" in t or "sars coronavirus" in t or m["organism"].lower() == "severe acute respiratory syndrome coronavirus": h = "human"
    else: h = "other"
    acc2host[m["accession"]] = h

def get(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "open-data-research/1.0"})
            with urllib.request.urlopen(req, timeout=90) as r:
                time.sleep(0.45); return r.read()
        except Exception as ex:
            print(" retry", ex, flush=True); time.sleep(2)
    raise RuntimeError(url)

uids = list(meta.keys())
spikes = {}
raw = open(f"{OUT}/raw/spike_cds_aa.fasta", "wb")
for i in range(0, len(uids), 100):
    chunk = uids[i:i+100]
    q = urllib.parse.urlencode({"db": "nuccore", "id": ",".join(chunk), "rettype": "fasta_cds_aa", "retmode": "text"})
    b = get(f"{E}/efetch.fcgi?{q}")
    raw.write(b)
    for rec in SeqIO.parse(io.StringIO(b.decode()), "fasta-pearson"):
        d = rec.description.lower()
        if not ("spike" in d or "surface glycoprotein" in d or "s glycoprotein" in d or " s protein" in d):
            continue
        if not (1100 <= len(rec.seq) <= 1400): continue
        prot = ""
        for k in ["protein_id=", "locus_tag="]:
            if k in rec.description: prot = rec.description.split(k)[1].split("]")[0]; break
        # link back to nucleotide accession via [location=..] or protein accession prefix
        spikes[rec.id] = (str(rec.seq), rec.description)
    print(f" batch {i}: spikes so far {len(spikes)}", flush=True)
raw.close()
print("raw spike CDS records:", len(spikes))
# map protein -> parent accession: efetch fasta_cds_aa headers carry [location=...] and protein id like YP_xxx; parent accession appears in description as [gbkey=CDS]... Use protein->parent map via elink? Instead: protein id prefix lookup in metadata title impossible.
# Simpler: the CDS fasta header contains [protein=spike glycoprotein] and protein id; parent = record order. We instead re-map via id table from epost history? Fallback: store by protein id, assign host by matching parent accession embedded in description if present.
import re
final = []
unmatched = 0
for pid, (seq, desc) in spikes.items():
    pm = re.match(r"lcl\|([A-Z0-9_.]+?)_prot_", pid)
    host = acc2host.get(pm.group(1), "other") if pm else "other"
    if not pm or host == "other": unmatched += 1
    final.append((pid, host, seq, desc))
print("unmatched parent:", unmatched)
# dedupe
seen, uniq = set(), []
for pid, h, seq, desc in final:
    if seq not in seen:
        seen.add(seq); uniq.append((pid, h, seq))
from collections import Counter
print("unique:", len(uniq), Counter(h for _, h, _ in uniq))
# SARS1 reference spike from genbank
gb = SeqIO.read(f"{OUT}/sars1_ref_NC_004718.gb", "genbank")
cands = [f for f in gb.features if f.type == "CDS" and ("surface" in str(f.qualifiers.get("product", [""])).lower() or "spike" in str(f.qualifiers.get("product", [""])).lower() or "S" in f.qualifiers.get("gene", []))]
s_cds = max(cands, key=lambda f: len(f))
ref_aa = s_cds.qualifiers["translation"][0]
with open("work/caseA_spikes.fasta", "w") as fh:
    fh.write(">SARS1_NC_004718_ref\n" + ref_aa + "\n")
    for pid, h, seq in uniq:
        fh.write(f">{pid}|{h}\n{seq}\n")
json.dump({pid: h for pid, h, _ in uniq}, open(f"{OUT}/caseA_spike_hosts.json", "w"), indent=0)
print("spikes fasta written:", len(uniq) + 1)
