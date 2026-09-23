# GATES_LOCKED — Project 13: Retrospective Outbreak Prediction + Pandemic Tools
Locked: 2026-09-23 (commit precedes ALL results commits; SHA256 of this file recorded in repo manifest)
Author: Builder 25 (for Udita Phookan)

## 0. Discipline (binding)
Open public data only. No wet lab. Free compute (2-core/1.9GB sandbox). No stubs, no
simulated data, no pseudocode-as-result, no placeholders. Real data, real numbers, real
figures. Honest verified/thin/missing counts in every report. Negative results are
preserved and reported; gates are never re-fit after outcomes are seen (no re-fishing).
A clearly-labelled EXPLORATORY section may follow pre-registered analysis but never
replaces it. The user authorized her connected accounts + internet tools for compute
ACCESS only; nothing is purchased and nothing is sent as her.

## 1. Core question
Using ONLY data publicly released before an outbreak, can per-site mutation-risk
predictions for a viral surface protein be scored — significantly better than chance and
than named prior-art baselines — against what actually evolved?

## 2. Case studies and locked date cutoffs
CASE A — SARS-CoV-2 (primary).
- PREDICTOR cutoff: only data publicly released on/before 2019-12-31.
  Predictor inputs: (i) Sarbecovirus genomes from NCBI GenBank with release date
  (PDAT) <= 2019-12-31; (ii) SARS-CoV-1 spike PDB structures released <= 2019-12-31
  (planned: 5XLR prefusion ectodomain; 2AJF RBD–ACE2 complex; release dates recorded
  in the data manifest). NO SARS-CoV-2 sequence, NO deep mutational scanning, nothing
  post-2019 enters any predictor.
- OUTCOME data (scoring only; never touches predictors): WHO VOC constellation
  defining spike substitutions — Alpha (B.1.1.7), Beta (B.1.351), Gamma (P.1),
  Delta (B.1.617.2), Omicron BA.1, BA.2 — from cov-lineages/constellations (public).
  Wuhan-Hu-1 / NC_045512 is used ONLY as the coordinate frame for scoring.

CASE B — Influenza A/H3N2 (transfer test).
- PREDICTOR cutoff: data released on/before 2009-12-31. Inputs: H3N2 HA sequences
  (NCBI, PDAT <= 2009-12-31); a pre-2010 HA crystal structure (release date recorded);
  antigenic epitope definitions as published pre-2010 (Wiley et al. 1981 sites A–E).
- OUTCOME data: Nextstrain seasonal-flu H3N2 HA tree; amino-acid substitutions that
  reached >50% global frequency with first appearance in 2010–2019.

## 3. Pre-registered method (predictor)
Per-site composite risk score, weights FIXED HERE (no fitting after outcomes seen):
  S(site) = mean of z-scored signals:
   s1 = amino-acid variability (Shannon entropy) across the pre-cutoff MSA
   s2 = structural exposure (Bio.PDB half-sphere exposure, inverted burial) from the
        pre-cutoff structure
   s3 = receptor-context proximity (min CA distance to receptor/contact footprint;
        SARS: ACE2 via 2AJF mapped onto spike; flu: sialic-acid receptor-binding site
        footprint as defined pre-2010)
   s4 = host-shift / adaptive-context indicator (Case A: site differs between
        bat-sarbecovirus consensus and human SARS-CoV-1 consensus; Case B: site lies
        in a pre-2010-defined antigenic epitope region A–E)
Sites lacking structural coverage are scored on available signals (mean of present z's)
and flagged; missing-data handling is part of the pre-registration.
PRE-REGISTERED PREDICTIONS = ranked site lists (Case A: top 30 RBD + top 20 NTD sites;
Case B: top 30 HA1 sites) written into PREREGISTRATION.md with method + scoring rule,
SHA256-hashed, and committed BEFORE any outcome data is downloaded or any scoring run.

## 4. Positive controls (must pass BEFORE any novel claim)
PC-A (flu): the pipeline, given only pre-2010 data, ranks >=5 of the 7 Koel et al.
  2013 antigenic-transition positions (H3 numbering 145,155,156,158,159,189,193) in the
  top quartile of scored HA1 sites. (Koel positions are outcome-side knowledge, used
  only as validation.)
PC-B (SARS): the pipeline's ACE2-contact annotation from 2AJF recovers >=80% of the
  published pre-2020 SARS-CoV-1 RBD–ACE2 contact-residue list.

## 5. Success gates (quantitative)
G1 = PC-A passes. G2 = PC-B passes.
G3 (primary, Case A): pre-registered predicted sites are enriched for VOC-defining
  spike substitutions vs a 10,000-draw label-permutation null over scored sites;
  one-sided p < 0.05.
G4 (primary, Case B): pre-registered HA1 risk sites are enriched for 2010–2019
  >50%-frequency substitutions vs the same permutation scheme; one-sided p < 0.05.
G5 (benchmark): site-level ROC-AUC of the composite score (positives = VOC-defining
  sites in RBD+NTD, Case A) exceeds the AUC of EVERY baseline below, with a bootstrap
  95% CI on the AUC difference vs the best baseline excluding 0:
   B1 = variability-only (s1); B2 = exposure-only (s2); B3 = Hie-style lightweight
   re-implementation (3-mer grammaticality log-likelihood under the pre-cutoff corpus
   + BLOSUM62 semantic-change distance); B4 = EVEscape-style independent-site baseline
   (pre-cutoff MSA amino-acid frequency x exposure).
Reported honestly either way: if G5 fails, the composite is reported as NOT beating
named prior-art baselines (a negative result, preserved).

## 6. Prior-art verdict
CROWDED — closest works named in PRIOR_ART.md (EVEscape / Thadani et al. Nature 2023;
Hie et al. Science 2021; Huddleston et al. eLife 2020; Obermeyer et al. Science 2022;
Maher et al. 2022). Differentiation: (i) fully pre-registered hindcast with git-hashed,
date-locked predictions (none of the named works pre-registered); (ii) transparent
laptop-scale composite (no GPU, no DMS, no pandemic-era data); (iii) dual-pathogen
transfer test; (iv) benchmarked vs lightweight re-implementations of named prior art;
(v) operational tooling built FROM the validated scores.

## 7. Pandemic tooling (built only from validated machinery)
T1 varscan-risk: per-site risk engine (MSA + optional PDB + receptor contact -> risk table)
T2 variant-interpret: annotate a variant's mutation list with risk tier, conservation,
   exposure, receptor contact, epitope membership
T3 target-impact: map mutations onto vaccine/antiviral target annotations (spike epitope
   classes; Mpro/RdRp functional sites) -> impact assessment
T4 sitrep: mutation-risk situation report generator (markdown) combining T2+T3

## 8. Seal criteria
Standalone gates commit (this file) first in git history; PREREGISTRATION.md committed
before outcome data acquisition; slice MANIFEST.sha256 covering data+results+code+paper
sources with byte-level sha256; fresh-sandbox restore + readback by orchestrator QC;
~20-page paper (PDF) following problem -> background/dataset research -> statistical
analysis -> results incl. honest negatives -> tools built from results; quantified
methodological contribution benchmarked against the named prior art above.
