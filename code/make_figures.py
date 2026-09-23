#!/usr/bin/env python3
"""All paper figures into results/figs/."""
import json, math, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Align import substitution_matrices
from collections import Counter
import os

os.makedirs("results/figs", exist_ok=True)
plt.rcParams.update({"figure.dpi": 150, "font.size": 9})

comp = json.load(open("data/predictor/composite_scores.json"))
compA = {int(k): v for k, v in comp["composite_A_sars1"].items()}
compB = {int(k): v for k, v in comp["composite_B_h3"].items()}
sigA = json.load(open("data/predictor/caseA_site_signals.json"))
sigB = json.load(open("data/predictor/caseB_site_signals.json"))
gates = json.load(open("results/gate_results.json"))
m12 = {int(k): int(v) for k, v in json.load(open("results/sars1_to_sars2_spike_map.json")).items()}
voc_pos = json.load(open("results/voc_spike_sites.json"))
all_voc = sorted(set(p for v in voc_pos.values() for p in v))

# ---- Fig 1: ROC curves ----
def zmap(d, invert=False):
    items = [(int(k), (-v if invert else v)) for k, v in d.items() if v is not None]
    xs = [v for _, v in items]; mu = sum(xs)/len(xs); sd = (sum((x-mu)**2 for x in xs)/len(xs))**0.5 or 1.0
    return {p: (v-mu)/sd for p, v in items}
z1 = zmap(sigA["s1_entropy"]); z2 = zmap(sigA["s2_hse_up"], invert=True)
msa = list(SeqIO.parse("work/caseA_spike_msa.fasta", "fasta"))
ref = next(r for r in msa if r.id == "SARS1_NC_004718_ref")
pos2col, pp = {}, 0
for j, c in enumerate(str(ref.seq)):
    if c != "-": pp += 1; pos2col[pp] = j
from Bio.Align import substitution_matrices as sm
bl = sm.load("BLOSUM62")
tri_counts = Counter()
for r in msa:
    sq = str(r.seq)
    for j in range(1, len(sq)-1):
        tri = sq[j-1:j+2]
        if "-" not in tri: tri_counts[tri] += 1
scored = [p for p in compA if 14 <= p <= 305 or 318 <= p <= 510]
voc_rbd = {p for p in all_voc if 319 <= p <= 541}
voc_ntd = {p for p in all_voc if 14 <= p <= 305}
voc_s2 = voc_rbd | voc_ntd
pos_s1 = {p for p, s2p in m12.items() if s2p in voc_s2}
sem, gram, consf = {}, {}, {}
for p in scored:
    j = pos2col[p]
    col = [str(r.seq)[j] for r in msa if str(r.seq)[j] != "-"]
    if not col: continue
    cons = Counter(col).most_common(1)[0][0]
    sem[p] = max([max(0.0, -bl[cons, a]) for a in set(col)-{cons}], default=0.0)
    tri = str(ref.seq)[j-1:j+2]
    gram[p] = math.log(tri_counts.get(tri, 0)+1) if "-" not in tri else 0.0
    consf[p] = 1.0 - Counter(col).most_common(1)[0][1]/len(col)
zsem, zgram, zc = zmap(sem), zmap(gram), zmap(consf)
models = {"PROMIS composite": {p: compA[p] for p in scored},
          "B1 variability": {p: z1.get(p, -9) for p in scored},
          "B2 exposure": {p: z2.get(p, -9) for p in scored},
          "B3 Hie-style": {p: zsem.get(p,0)+zgram.get(p,0) for p in scored},
          "B4 EVEscape-style": {p: zc.get(p,0)+z2.get(p,0) for p in scored}}
y = np.array([1 if p in pos_s1 else 0 for p in scored])
fig, ax = plt.subplots(figsize=(5.2, 4.4))
aucs = {}
for name, m in models.items():
    s = np.array([m[p] for p in scored])
    order = np.argsort(-s)
    ys = y[order]
    tpr = np.cumsum(ys)/ys.sum(); fpr = np.cumsum(1-ys)/(len(ys)-ys.sum())
    auc_v = np.trapz(np.concatenate([[0], tpr]), np.concatenate([[0], fpr]))
    aucs[name] = auc_v
    ax.plot(np.concatenate([[0], fpr]), np.concatenate([[0], tpr]),
            label=f"{name} (AUC={auc_v:.3f})", lw=2 if "PROMIS" in name else 1)
ax.plot([0,1],[0,1],"k--",lw=0.7)
ax.set_xlabel("FPR"); ax.set_ylabel("TPR")
ax.set_title("Site-level VOC-mutation prediction, SARS-CoV-2 RBD+NTD\n(predictors: pre-2020 data only)")
ax.legend(fontsize=7.5, loc="lower right")
fig.tight_layout(); fig.savefig("results/figs/fig1_roc_caseA.png"); plt.close(fig)
print("fig1 aucs:", {k: round(v,3) for k,v in aucs.items()})

# ---- Fig 2: permutation nulls G3/G4 ----
fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6))
# regenerate nulls deterministically (same code as scoring)
rng = np.random.default_rng(13)
pred_rbd_s2 = {m12[p] for p in [int(x) for x in re.search(r"Top-30 RBD \(318-510\): ([0-9, ]+)", open("PREREGISTRATION.md").read()).group(1).split(",")] if p in m12}
pred_ntd_s2 = {m12[p] for p in [int(x) for x in re.search(r"Top-20 NTD \(14-305\): ([0-9, ]+)", open("PREREGISTRATION.md").read()).group(1).split(",")] if p in m12}
scored_rbd_s2 = [m12[p] for p in scored if 318 <= p <= 510 and p in m12]
scored_ntd_s2 = [m12[p] for p in scored if 14 <= p <= 305 and p in m12]
nullA = []
for _ in range(10000):
    r = rng.choice(scored_rbd_s2, size=len(pred_rbd_s2), replace=False)
    n = rng.choice(scored_ntd_s2, size=len(pred_ntd_s2), replace=False)
    nullA.append(len(set(r) & voc_rbd) + len(set(n) & voc_ntd))
obsA = gates["G3"]["observed_hits"]
axes[0].hist(nullA, bins=range(0, max(nullA)+2), color="#9ecae1", edgecolor="white")
axes[0].axvline(obsA, color="red", lw=2, label=f"observed = {obsA}")
axes[0].set_title(f"G3 null (10,000 permutations)\np = {gates['G3']['p']:.4f}")
axes[0].set_xlabel("VOC sites hit by random site set"); axes[0].legend()
est = set(gates["G4"]["established"])
pred_ha1 = set(gates["G4"]["pred"])
scored_ha1 = list(compB.keys())
nullB = [len(set(rng.choice(scored_ha1, size=len(pred_ha1), replace=False)) & est) for _ in range(10000)]
obsB = gates["G4"]["observed_hits"]
axes[1].hist(nullB, bins=range(0, max(nullB)+2), color="#a1d99b", edgecolor="white")
axes[1].axvline(obsB, color="red", lw=2, label=f"observed = {obsB}")
axes[1].set_title(f"G4 null (10,000 permutations)\np = {gates['G4']['p']:.2e}")
axes[1].set_xlabel("established-substitution sites hit by random set"); axes[1].legend()
fig.tight_layout(); fig.savefig("results/figs/fig2_permutations.png"); plt.close(fig)
print("fig2 done")

# ---- Fig 3: RBD score profile + VOC markers ----
fig, ax = plt.subplots(figsize=(9.5, 3.4))
xs = sorted(p for p in compA if 318 <= p <= 510)
ax.plot(xs, [compA[p] for p in xs], color="#3182bd", lw=1.2, label="PROMIS score (pre-2020 data)")
hit_s1 = {p for p, s in m12.items() if s in voc_rbd}
for p in xs:
    if p in hit_s1:
        ax.axvline(p, color="red", alpha=0.55, lw=1)
top30 = set(int(x) for x in re.search(r"Top-30 RBD \(318-510\): ([0-9, ]+)", open("PREREGISTRATION.md").read()).group(1).split(","))
ax.scatter([p for p in xs if p in top30], [compA[p] for p in xs if p in top30],
           s=18, facecolor="none", edgecolor="black", label="pre-registered top-30")
ax.plot([], [], color="red", lw=1, label="actual VOC-defining site (2020-22)")
ax.set_xlabel("spike position (SARS-CoV-1 numbering)"); ax.set_ylabel("composite z")
ax.set_title("Case A: predicted RBD risk profile vs actual VOC mutations")
ax.legend(fontsize=7.5); fig.tight_layout(); fig.savefig("results/figs/fig3_rbd_profile.png"); plt.close(fig)
print("fig3 done")

# ---- Fig 4: HA1 profile with Koel + established ----
fig, ax = plt.subplots(figsize=(9.5, 3.4))
xs = sorted(compB)
ax.plot(xs, [compB[p] for p in xs], color="#31a354", lw=1.2, label="PROMIS score (pre-2010 data)")
KOEL = [145,155,156,158,159,189,193]
for p in KOEL: ax.axvline(p, color="purple", alpha=0.5, lw=1)
for p in est: ax.axvline(p, color="red", alpha=0.25, lw=1)
ax.plot([], [], color="purple", lw=1, label="Koel antigenic positions (positive control)")
ax.plot([], [], color="red", lw=1, alpha=0.4, label="substitutions established 2010s (Nextstrain)")
top30B = set(gates["G4"]["pred"])
ax.scatter([p for p in xs if p in top30B], [compB[p] for p in xs if p in top30B],
           s=18, facecolor="none", edgecolor="black", label="pre-registered top-30")
ax.set_xlabel("HA1 position (H3 mature numbering)"); ax.set_ylabel("composite z")
ax.set_title("Case B: predicted H3N2 HA1 risk profile vs actual post-2010 evolution")
ax.legend(fontsize=7.5); fig.tight_layout(); fig.savefig("results/figs/fig4_ha1_profile.png"); plt.close(fig)
print("fig4 done")

# ---- Fig 5: signal heatmap top RBD sites ----
z3 = zmap(sigA["s3_neg_ace2_dist"]); z4 = zmap(sigA["s4_hostshift"])
top30_list = [int(x) for x in re.search(r"Top-30 RBD \(318-510\): ([0-9, ]+)", open("PREREGISTRATION.md").read()).group(1).split(",")]
M = np.array([[z1.get(p, np.nan), z2.get(p, np.nan), z3.get(p, np.nan), z4.get(p, np.nan), compA[p]] for p in top30_list])
fig, ax = plt.subplots(figsize=(6.2, 5.6))
im = ax.imshow(M, aspect="auto", cmap="viridis")
ax.set_yticks(range(len(top30_list))); ax.set_yticklabels([str(p) for p in top30_list], fontsize=7)
ax.set_xticks(range(5)); ax.set_xticklabels(["variability", "exposure", "ACE2 prox.", "host shift", "PROMIS"], rotation=30, ha="right")
ax.set_title("Signal decomposition, pre-registered top-30 RBD sites (z-scores)")
ax.set_ylabel("spike position (SARS-CoV-1)")
fig.colorbar(im, ax=ax, shrink=0.8)
fig.tight_layout(); fig.savefig("results/figs/fig5_signal_heatmap.png"); plt.close(fig)
print("fig5 done")

# ---- Fig 6: flu trajectories of established predicted sites ----
tree = json.load(open("data/outcome/nextstrain_h3n2_ha_12y.json"))["tree"]
fr_raw = json.load(open("data/outcome/nextstrain_h3n2_tip_frequencies.json"))
pivots = fr_raw["pivots"]
strain_freqs = {k: v["frequencies"] for k, v in fr_raw.items() if isinstance(v, dict) and "frequencies" in v}
def walk(n, acc):
    acc.append(n)
    for c in n.get("children", []): walk(c, acc)
nodes = []; walk(tree, nodes)
def tips(n):
    if not n.get("children"):
        nm = n.get("name", "")
        return [nm] if nm in strain_freqs else []
    t = []
    for c in n["children"]: t += tips(c)
    return t
trajs = {}
for n in nodes:
    d = n.get("node_attrs", {}).get("num_date", {}).get("value")
    muts = (n.get("branch_attrs", {}).get("mutations", {}) or {}).get("HA1", [])
    if not muts or d is None or not (2010 <= d < 2020): continue
    for mut in muts:
        mm = re.match(r"([A-Z])(\d+)([A-Z])", mut)
        if not mm: continue
        pos = int(mm.group(2))
        if pos not in (set(gates["G4"]["hits"])): continue
        tp = tips(n)
        if not tp: continue
        fr = np.zeros(len(pivots))
        for t in tp:
            v = strain_freqs.get(t)
            if v: fr += np.array(v[:len(pivots)] + [0.0]*max(0, len(pivots)-len(v)))
        trajs[f"{mut} (d~{d:.1f})"] = fr
fig, ax = plt.subplots(figsize=(8.6, 4.2))
for label, fr in sorted(trajs.items(), key=lambda kv: -kv[1].max())[:8]:
    ax.plot(pivots, fr, lw=1.3, label=label)
ax.set_xlabel("year"); ax.set_ylabel("estimated global frequency")
ax.set_title("Case B: pre-registered sites whose substitutions swept to high frequency after 2010")
ax.legend(fontsize=7, ncol=2)
fig.tight_layout(); fig.savefig("results/figs/fig6_flu_trajectories.png"); plt.close(fig)
print("fig6 done:", len(trajs), "trajectories")
