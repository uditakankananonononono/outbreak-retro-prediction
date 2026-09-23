# outbreak-retro-prediction
Pre-registered outbreak mutation/site predictions using only pre-outbreak data, scored
against actual evolution, plus pandemic-response tooling built from the validated scores.
Project 13 of the 13-project open-data computational-biology portfolio.

## Headline result
Frozen, hash-locked site predictions (PREREGISTRATION.md, sha256 d33a3f28..., committed
BEFORE outcome data) were enriched for realized evolution on two pathogens:
- SARS-CoV-2: 7 VOC-defining spike sites hit (E484, N501, Q493, G496, L452, D405, T213)
  vs 2.1 expected by permutation, p = 0.0016 (G3 PASS)
- Influenza H3N2: 14 of 33 sites whose substitutions swept 2010-2019 vs 3.0 expected,
  p < 0.0001 (G4 PASS); positive controls 7/7 Koel (PC-A) and 87.5% ACE2 contacts (PC-B)
- Honest negative: composite AUC 0.722 does NOT beat EVEscape-style baseline 0.724 (G5 FAIL, preserved)

## Repo map
- GATES_LOCKED.md / PRIOR_ART.md / PREREGISTRATION.md - discipline, prior art, frozen predictions
- code/ - pipeline in run order: fetch_predictor_data.py -> extract_spikes_v2.py ->
  build_caseA_scores.py / build_caseB_scores.py -> make_prereg.py -> fetch_outcome_data.py ->
  score_gates.py -> make_figures.py -> make_paper.py -> make_manifest.py
- data/ - predictor/ + outcome/ payloads with sha256 manifests (raw API responses in raw/)
- results/ - gate_results.json, figures/, voc sites, coordinate map, BA.1 demo sitrep
- tools/ - T1 varscan_risk.py, T2 variant_interpret.py, T3 target_impact.py, T4 sitrep.py
- paper/outbreak_retro_prediction.pdf - 17-page manuscript
- MANIFEST.sha256 - slice byte-lock (81 files) for fresh-sandbox restore QC

## Quickstart
python3 tools/variant_interpret.py --case sarscov2 --mutations S:E484K,S:N501Y,S:Q493R
python3 tools/variant_interpret.py --case h3n2 --mutations HA1:N145D,HA1:K160T
python3 tools/target_impact.py --mutations S:E484K,S:P681H
python3 tools/sitrep.py --name "NewVariant" --mutations S:E484K,S:N501Y --out report.md
Full reproduction: run code/ scripts in the order above (NCBI eutils, RCSB, Nextstrain;
anonymous public endpoints only; muscle 5.3 binary required for alignments).
