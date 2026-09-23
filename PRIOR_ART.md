# PRIOR ART SWEEP — Project 13 (run 2026-09-23, before gates locked)
Verdict: CROWDED (build, with explicit differentiation — not DONE, not CLEAR)

Closest works:
1. Thadani, Gurev et al. "Learning from prepandemic data to forecast viral escape."
   Nature 622, 2023. https://www.nature.com/articles/s41586-023-06617-0
   EVEscape: deep generative model trained on pre-pandemic coronavirus sequences +
   DMS + structure; forecasts SARS-CoV-2 escape. CLOSEST. Differentiators: they did not
   pre-register site lists; model needs heavy training; scores escape, not realized
   VOC emergence; no transfer test to a second pathogen.
2. Hie, Zhong, Berger, Bryson. "Learning the language of viral evolution and escape."
   Science 371, 2021. https://www.science.org/doi/10.1126/science.abd7331
   Language models (grammaticality + semantic change) trained on pre-pandemic data
   predict escape. We re-implement a lightweight version (baseline B3) and benchmark.
3. Huddleston et al. "Integrating genotypes and phenotypes improves long-term forecasts
   of seasonal influenza A/H3N2 evolution." eLife 9:e60067, 2020.
   https://elifesciences.org/articles/60067 — flu fitness forecasting from within-season
   data (different regime: they forecast near-term using circulating diversity; we
   hindcast from a hard pre-outbreak cutoff).
4. Obermeyer et al. "Analysis of 6.4 million SARS-CoV-2 genomes identifies mutations
   implicated in fitness." Science 376, 2022 — during-pandemic data, not pre-outbreak.
5. Maher et al. "Predicting the mutational drivers of future SARS-CoV-2 variants of
   concern." Science Advances 8, 2022 — PyR0 fitness model on pandemic genomes.
6. Neher, Bedford. "Nextflu: real-time tracking of seasonal influenza A/H3N2 evolution"
   / LBI (eLife 2015 era) — tree-based fitness; needs circulating-diversity trees.
7. Deep-learning SARS-CoV-2 evolution predictors, 2023–2026 (e.g. Nature Signal Transduct
   Target Ther 2024 s41392-024-02066-x; Nature Microbiol 2026 s41564-026-02377-5) —
   all train on pandemic-era data; not pre-outbreak-only.

Our claimed methodological contribution (quantified in results): a fully pre-registered,
interpretable, laptop-scale composite site-risk score (PROMIS: Pre-outbreak Mutation
Risk Score) benchmarked site-by-site against lightweight re-implementations of named
prior art (B1–B4 in GATES_LOCKED.md) on two pathogens with hard date cutoffs, with
permutation statistics — plus operational tooling (T1–T4) built from the validated score.
