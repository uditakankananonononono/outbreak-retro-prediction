#!/usr/bin/env python3
"""Fetch PREDICTOR-side data under locked date cutoffs. Writes data/predictor/ + manifest.
Case A cutoff: 2019-12-31 (sarbecovirus genomes, SARS-CoV-1 structures).
Case B cutoff: 2009-12-31 (H3N2 HA sequences, pre-2010 HA structure).
All sources anonymous public APIs. Sleeps respect NCBI 3 req/s."""
import json, subprocess, time, hashlib, os, sys, urllib.parse, urllib.request

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = "data/predictor"
os.makedirs(f"{OUT}/raw", exist_ok=True)

def get(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"open-data-research/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                time.sleep(0.45)
                return r.read()
        except Exception as ex:
            print(f"  retry {i+1} {ex}", flush=True); time.sleep(2+2*i)
    raise RuntimeError("failed: "+url)

def esearch(term, retmax=10000):
    q = urllib.parse.urlencode({"db":"nuccore","term":term,"retmode":"json","retmax":retmax})
    d = json.loads(get(f"{E}/esearch.fcgi?{q}"))
    return d["esearchresult"]["idlist"]

def esummary(ids):
    out = {}
    for i in range(0, len(ids), 200):
        chunk = ids[i:i+200]
        q = urllib.parse.urlencode({"db":"nuccore","id":",".join(chunk),"retmode":"json","version":"2.0"})
        d = json.loads(get(f"{E}/esummary.fcgi?{q}"))
        with open(f"{OUT}/raw/esummary_{i}.json","w") as f: json.dump(d,f)
        for uid in d["result"]["uids"]:
            r = d["result"][uid]
            out[uid] = {"accession": r.get("accessionversion",""), "title": r.get("title",""),
                        "organism": r.get("organism",""), "createdate": r.get("createdate",""),
                        "updatedate": r.get("updatedate",""), "slen": r.get("slen",0)}
        print(f"  esummary {i+len(chunk)}/{len(ids)}", flush=True)
    return out

def efetch_fasta(ids, path):
    n = 0
    with open(path, "wb") as fh:
        for i in range(0, len(ids), 150):
            chunk = ids[i:i+150]
            q = urllib.parse.urlencode({"db":"nuccore","id":",".join(chunk),"rettype":"fasta","retmode":"text"})
            b = get(f"{E}/efetch.fcgi?{q}")
            fh.write(b); n += b.count(b">")
            print(f"  efetch {n} seqs", flush=True)
    return n

# ---- Case A: sarbecovirus genomes, PDAT <= 2019-12-31 ----
print("== Case A: sarbecovirus esearch ==", flush=True)
idsA = esearch('Sarbecovirus[Organism] AND 20000:32000[SLEN] AND ("1900/01/01"[PDAT] : "2019/12/31"[PDAT])')
json.dump(idsA, open(f"{OUT}/raw/esearch_sarbecovirus_ids.json","w"))
print(f"  {len(idsA)} ids", flush=True)
metaA = esummary(idsA)
# defensible filter: GenBank CreateDate (first public release) <= 2019-12-31
keepA = {u:m for u,m in metaA.items() if m["createdate"] and m["createdate"] <= "2019/12/31"}
dropA = {u:m for u,m in metaA.items() if u not in keepA}
print(f"  kept {len(keepA)} / dropped {len(dropA)} by CreateDate", flush=True)
for u,m in list(dropA.items())[:10]: print("   DROP", m["accession"], m["createdate"], m["organism"][:50], flush=True)
nA = efetch_fasta(list(keepA.keys()), f"{OUT}/sarbecovirus_genomes.fasta")
json.dump(keepA, open(f"{OUT}/sarbecovirus_metadata.json","w"), indent=1)
print(f"Case A: {nA} genomes fetched", flush=True)

# ---- Case B: H3N2 HA, PDAT <= 2009-12-31 ----
print("== Case B: H3N2 HA esearch ==", flush=True)
idsB = esearch('Influenza A virus[Organism] AND H3N2[All Fields] AND 1650:1800[SLEN] AND ("1968/01/01"[PDAT] : "2009/12/31"[PDAT]) AND (hemagglutinin[All Fields] OR "segment 4"[All Fields])')
json.dump(idsB, open(f"{OUT}/raw/esearch_h3n2_ids.json","w"))
print(f"  {len(idsB)} ids", flush=True)
metaB = esummary(idsB)
keepB = {u:m for u,m in metaB.items() if m["createdate"] and m["createdate"] <= "2009/12/31"}
print(f"  kept {len(keepB)} / dropped {len(metaB)-len(keepB)} by CreateDate", flush=True)
nB = efetch_fasta(list(keepB.keys()), f"{OUT}/h3n2_ha.fasta")
json.dump(keepB, open(f"{OUT}/h3n2_metadata.json","w"), indent=1)
print(f"Case B: {nB} HA sequences fetched", flush=True)

# ---- SARS-CoV-1 reference (2003, pre-cutoff): genbank + fasta ----
gb = get(f"{E}/efetch.fcgi?db=nuccore&id=NC_004718.3&rettype=gbwithparts&retmode=text")
open(f"{OUT}/sars1_ref_NC_004718.gb","wb").write(gb)
fa = get(f"{E}/efetch.fcgi?db=nuccore&id=NC_004718.3&rettype=fasta&retmode=text")
open(f"{OUT}/sars1_ref_NC_004718.fasta","wb").write(fa)

# ---- PDB entries: release dates + mmCIF ----
for pdb in ["5XLR","2AJF"]:
    info = json.loads(get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb}"))
    rel = info.get("rcsb_accession_info",{}).get("initial_release_date","?")
    title = info.get("struct",{}).get("title","")
    print(f"PDB {pdb}: released {rel} | {title[:80]}", flush=True)
    cif = get(f"https://files.rcsb.org/download/{pdb}.cif")
    open(f"{OUT}/{pdb}.cif","wb").write(cif)
    open(f"{OUT}/raw/{pdb}_entry.json","w").write(json.dumps({"release":rel,"title":title}))

# ---- checksums ----
man = {}
for root,_,files in os.walk(OUT):
    for fn in files:
        p = os.path.join(root,fn)
        man[p] = hashlib.sha256(open(p,'rb').read()).hexdigest()
json.dump(man, open("data/predictor_manifest.json","w"), indent=1)
print(f"manifest: {len(man)} files byte-locked", flush=True)
