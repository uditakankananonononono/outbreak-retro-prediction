#!/usr/bin/env python3
"""Score locked gates G1-G5 against outcome data. Runs ONLY after PREREGISTRATION commit."""
import json, math, re, subprocess
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
import numpy as np

rng = np.random.default_rng(13)
OUTJ = {}

# ---------- shared data ----------
comp = json.load(open("data/predictor/composite_scores.json"))
compA = {int(k): v for k, v in comp["composite_A_sars1"].items()}
compB = {int(k): v for k, v in comp["composite_B_h3"].items()}
sigA = json.load(open("data/predictor/caseA_site_signals.json"))
sigB = json.load(open("data/predictor/caseB_site_signals.json"))
prereg = open("PREREGISTRATION.md").read()
top_rbd = [int(x) for x in re.search(r"Top-30 RBD \(318-510\): ([0-9, ]+)", prereg).group(1).split(",")]
top_ntd = [int(x) for x in re.search(r"Top-20 NTD \(14-305\): ([0-9, ]+)", prereg).group(1).split(",")]
top_ha1 = [int(x) for x in re.search(r"Top-30 HA1: ([0-9, ]+)", prereg).group(1).split(",")]

# ---------- SARS1 -> SARS2 spike coordinate map (SARS2 ref = coordinate frame only) ----------
gb2 = SeqIO.read("data/outcome/sars2_ref_NC_045512.gb", "genbank")
s2_cds = next(f for f in gb2.features if f.type == "CDS" and "surface" in str(f.qualifiers.get("product", [""])).lower())
sars2_spike = s2_cds.qualifiers["translation"][0]
sars1_spike = str(SeqIO.read("data/predictor/sars1_ref_NC_004718.fasta", "fasta").seq)[21491:25259]
sars1_spike_aa = str(__import__("Bio.Seq", fromlist=["Seq"]).Seq(sars1_spike).translate())
al = Align.PairwiseAligner()
al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.open_gap_score = -10; al.extend_gap_score = -0.5; al.mode = "global"
a = al.align(sars1_spike_aa, sars2_spike)[0]
m12 = {}
for (s1, e1), (s2, e2) in zip(a.aligned[0], a.aligned[1]):
    for k in range(e1 - s1):
        m12[int(s1 + k + 1)] = int(s2 + k + 1)  # SARS1 pos -> SARS2 pos
json.dump(m12, open("results/sars1_to_sars2_spike_map.json", "w"))

# ---------- VOC sites from constellations (SARS2 numbering) ----------
voc_pos = {}
for voc in ["Alpha", "Beta", "Gamma", "Delta", "Omicron_BA1", "Omicron_BA2"]:
    d = json.load(open(f"data/outcome/constellation_{voc}.json"))
    ps = set()
    for site in d.get("sites", []):
        m = re.match(r"(?:spike|s):([A-Z*])?(\d+)", str(site), re.I)
        if m: ps.add(int(m.group(2)))
    voc_pos[voc] = sorted(ps)
    print(voc, len(ps), "spike sites")
all_voc = sorted(set().union(*[set(v) for v in voc_pos.values()]))
json.dump(voc_pos, open("results/voc_spike_sites.json", "w"), indent=1)

# ---------- PC-B: ACE2 contact recovery (Li et al. 2005 Science 309:1864 list) ----------
LI2005 = [402, 424, 426, 436, 440, 442, 472, 473, 475, 479, 482, 484, 486, 487, 488, 491]
s3 = {int(k): v for k, v in sigA["s3_neg_ace2_dist"].items()}
contacts = {p for p, v in s3.items() if v > -8.0}
pcb_hits = [p for p in LI2005 if p in contacts]
OUTJ["PC_B"] = {"published_n": len(LI2005), "recovered": len(pcb_hits),
                "missing": [p for p in LI2005 if p not in contacts],
                "pass": len(pcb_hits) >= math.ceil(0.8 * len(LI2005))}
print("PC-B:", OUTJ["PC_B"])

# ---------- PC-A / G1: Koel positions in top quartile (Case B) ----------
KOEL = [145, 155, 156, 158, 159, 189, 193]
b_scores = np.array(sorted(compB.values()))
q75 = np.quantile(b_scores, 0.75)
koel_hits = [p for p in KOEL if compB.get(p, -9) >= q75]
OUTJ["PC_A"] = {"koel": KOEL, "in_top_quartile": koel_hits, "pass": len(koel_hits) >= 5}
print("PC-A:", OUTJ["PC_A"])

# ---------- G3: Case A enrichment (permutation) ----------
scored_rbd = [p for p in compA if 318 <= p <= 510]
scored_ntd = [p for p in compA if 14 <= p <= 305]
pred_rbd_s2 = {m12[p] for p in top_rbd if p in m12}
pred_ntd_s2 = {m12[p] for p in top_ntd if p in m12}
voc_rbd = {p for p in all_voc if 319 <= p <= 541}
voc_ntd = {p for p in all_voc if 14 <= p <= 305}
obs = len(pred_rbd_s2 & voc_rbd) + len(pred_ntd_s2 & voc_ntd)
scored_rbd_s2 = [m12[p] for p in scored_rbd if p in m12]
scored_ntd_s2 = [m12[p] for p in scored_ntd if p in m12]
null = []
for _ in range(10000):
    r = rng.choice(scored_rbd_s2, size=len(pred_rbd_s2), replace=False)
    n = rng.choice(scored_ntd_s2, size=len(pred_ntd_s2), replace=False)
    null.append(len(set(r) & voc_rbd) + len(set(n) & voc_ntd))
null = np.array(null)
p_g3 = (1 + (null >= obs).sum()) / 10001
OUTJ["G3"] = {"pred_rbd_n": len(pred_rbd_s2), "pred_ntd_n": len(pred_ntd_s2),
              "voc_rbd": sorted(voc_rbd), "voc_ntd": sorted(voc_ntd),
              "hits_rbd": sorted(pred_rbd_s2 & voc_rbd), "hits_ntd": sorted(pred_ntd_s2 & voc_ntd),
              "observed_hits": obs, "null_mean": float(null.mean()), "p": float(p_g3),
              "pass": bool(p_g3 < 0.05)}
print("G3:", OUTJ["G3"]["observed_hits"], "hits, p =", p_g3)

# ---------- G4: Case B enrichment vs Nextstrain-established HA1 subs 2010-2019 ----------
tree = json.load(open("data/outcome/nextstrain_h3n2_ha_12y.json"))["tree"]
freqs_raw = json.load(open("data/outcome/nextstrain_h3n2_tip_frequencies.json"))
pivots = freqs_raw.get("pivots") or freqs_raw.get("data", {}).get("pivots")
strain_freqs = {}
for k, v in freqs_raw.items():
    if isinstance(v, dict) and "frequencies" in v:
        strain_freqs[k] = v["frequencies"]
NP = len(pivots)
def freqvec(t):
    a = strain_freqs.get(t)
    if a is None: return None
    v = np.array(a[:NP], dtype=float)
    if len(v) < NP: v = np.pad(v, (0, NP - len(v)))
    return v
print("freq strains:", len(strain_freqs), "pivots:", len(pivots) if pivots else None)

nodes = []
def walk(n, parent=None):
    nodes.append(n); n["_parent"] = parent
    for c in n.get("children", []): walk(c, n)
walk(tree)
for n in nodes:
    n["_tips"] = []
def collect(n):
    if not n.get("children"):
        nm = n.get("name") or n.get("node_attrs", {}).get("strain", "")
        n["_tips"] = [nm] if nm in strain_freqs else []
        return n["_tips"]
    t = []
    for c in n["children"]: t += collect(c)
    n["_tips"] = t; return t
collect(tree)
pivot_idx = [i for i, pv in enumerate(pivots) if 2010.0 <= pv < 2020.0]
est = {}  # pos -> max freq
for n in nodes:
    date = n.get("node_attrs", {}).get("num_date", {}).get("value")
    muts = (n.get("branch_attrs", {}).get("mutations", {}) or {}).get("HA1", [])
    if not muts or date is None or not (2010.0 <= date < 2020.0): continue
    if not n["_tips"]: continue
    fr = np.zeros(NP)
    for t in n["_tips"]:
        v = freqvec(t)
        if v is not None: fr += v
    mx = fr.max() if NP else 0.0
    for mut in muts:
        m = re.match(r"[A-Z](\d+)[A-Z]", mut)
        if m:
            p = int(m.group(1))
            est[p] = max(est.get(p, 0.0), float(mx))
established = {p for p, f in est.items() if f >= 0.5}
print("established HA1 positions 2010-19 (>50%):", sorted(established))
scored_ha1 = list(compB.keys())
pred_set = set(top_ha1)
obsB = len(pred_set & established)
nullB = []
for _ in range(10000):
    r = rng.choice(scored_ha1, size=len(pred_set), replace=False)
    nullB.append(len(set(r) & established))
nullB = np.array(nullB)
p_g4 = (1 + (nullB >= obsB).sum()) / 10001
OUTJ["G4"] = {"established": sorted(established), "pred": sorted(pred_set),
              "hits": sorted(pred_set & established), "observed_hits": obsB,
              "null_mean": float(nullB.mean()), "p": float(p_g4), "pass": bool(p_g4 < 0.05)}
print("G4:", obsB, "hits, p =", p_g4)

# ---------- G5: AUC vs baselines (Case A, RBD+NTD) ----------
def zmap(d, invert=False):
    items = [(int(k), (-v if invert else v)) for k, v in d.items() if v is not None]
    xs = [v for _, v in items]; mu = sum(xs) / len(xs)
    sd = (sum((x - mu) ** 2 for x in xs) / len(xs)) ** 0.5 or 1.0
    return {p: (v - mu) / sd for p, v in items}

scored = [p for p in compA if 14 <= p <= 305 or 318 <= p <= 510]
pos_s1 = {p for p, s2p in m12.items() if s2p in (voc_rbd | voc_ntd)}
y = np.array([1 if p in pos_s1 else 0 for p in scored])

def auc(scores):
    order = np.argsort(scores)
    ranks = np.empty(len(scores)); ranks[order] = np.arange(1, len(scores) + 1)
    # average ties
    s_sorted = scores[order]
    i = 0
    while i < len(s_sorted):
        j = i
        while j + 1 < len(s_sorted) and s_sorted[j + 1] == s_sorted[i]: j += 1
        if j > i:
            avg = (ranks[order[i]] + ranks[order[j]]) / 2
            ranks[order[i:j + 1]] = avg
        i = j + 1
    n1 = y.sum(); n0 = len(y) - n1
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))

z1 = zmap(sigA["s1_entropy"]); z2 = zmap(sigA["s2_hse_up"], invert=True)
# B3: Hie-style: semantic-change potential (max -BLOSUM sim of observed alts vs consensus) + 3-mer context richness
msa = list(SeqIO.parse("work/caseA_spike_msa.fasta", "fasta"))
ref = next(r for r in msa if r.id == "SARS1_NC_004718_ref")
pos2col, pp = {}, 0
for j, c in enumerate(str(ref.seq)):
    if c != "-": pp += 1; pos2col[pp] = j
bl = substitution_matrices.load("BLOSUM62")
from collections import Counter
sem, gram = {}, {}
tri_counts = Counter()
for r in msa:
    sq = str(r.seq)
    for j in range(1, len(sq) - 1):
        tri = sq[j - 1:j + 2]
        if "-" not in tri: tri_counts[tri] += 1
for p in scored:
    j = pos2col[p]
    col = [str(r.seq)[j] for r in msa if str(r.seq)[j] != "-"]
    if not col: continue
    cons = Counter(col).most_common(1)[0][0]
    alts = set(col) - {cons}
    sem[p] = max([max(0.0, -bl[cons, a]) for a in alts], default=0.0)
    tri = str(ref.seq)[j - 1:j + 2]
    gram[p] = math.log(tri_counts.get(tri, 0) + 1) if "-" not in tri else 0.0
zsem = zmap(sem); zgram = zmap(gram)
B3 = {p: zsem.get(p, 0) + zgram.get(p, 0) for p in scored}
# B4: EVEscape-style independent-site: (1 - consensus freq) + exposure
consf = {}
for p in scored:
    j = pos2col[p]
    col = [str(r.seq)[j] for r in msa if str(r.seq)[j] != "-"]
    if col: consf[p] = 1.0 - Counter(col).most_common(1)[0][1] / len(col)
zc = zmap(consf)
B4 = {p: zc.get(p, 0) + z2.get(p, 0) for p in scored}

models = {"COMPOSITE": {p: compA[p] for p in scored},
          "B1_variability": {p: z1.get(p, -9) for p in scored},
          "B2_exposure": {p: z2.get(p, -9) for p in scored},
          "B3_hie_style": B3, "B4_evescape_style": B4}
aucs = {name: auc(np.array([m[p] for p in scored])) for name, m in models.items()}
print("AUCs:", aucs)
best_base = max((v, k) for k, v in aucs.items() if k != "COMPOSITE")[1]
diff_obs = aucs["COMPOSITE"] - aucs[best_base]
boot = []
idx = np.arange(len(scored))
for _ in range(10000):
    b = rng.choice(idx, size=len(idx), replace=True)
    if y[b].sum() in (0, len(b)): continue
    ac = auc(np.array([models["COMPOSITE"][scored[i]] for i in b]))
    ab = auc(np.array([models[best_base][scored[i]] for i in b]))
    boot.append(ac - ab)
boot = np.array(boot)
ci = np.quantile(boot, [0.025, 0.975])
OUTJ["G5"] = {"aucs": aucs, "best_baseline": best_base, "auc_diff": float(diff_obs),
              "ci95": [float(ci[0]), float(ci[1])], "pass": bool(ci[0] > 0)}
print("G5:", OUTJ["G5"])
OUTJ["VOC_sites"] = voc_pos
json.dump(OUTJ, open("results/gate_results.json", "w"), indent=1)
print("gate results written")
