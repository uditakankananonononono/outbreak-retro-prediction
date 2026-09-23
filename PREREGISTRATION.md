# PREREGISTRATION - Project 13 retrospective outbreak prediction
Generated: 2026-09-23T15:54:06.346410Z. Committed BEFORE any outcome data download or scoring.

## Method (as locked in GATES_LOCKED.md)
Composite site-risk score PROMIS = unweighted mean of z-scored signals s1 variability,
s2 exposure (HSE inverted), s3 receptor proximity (negative min-CA distance), s4 context
(Case A bat-vs-human consensus shift; Case B pre-2010 epitope membership). Site eligible if
MSA coverage >= 0.5 and >= 3 of 4 signals present. No weights fitted anywhere.

## Scoring rules (locked)
- G3: Case A enrichment of VOC-defining spike substitutions (WHO VOCs Alpha, Beta, Gamma,
  Delta, Omicron BA.1+BA.2, from cov-lineages constellations) among predicted sites vs
  10,000-draw label-permutation null over all scored RBD+NTD sites; one-sided p<0.05.
- G4: Case B enrichment of HA1 substitutions reaching >50% global frequency (first
  appearance 2010-2019, Nextstrain 12y tree + tip frequencies) among predicted sites; same null.
- G5: site-level ROC-AUC of composite (Case A, positives = VOC sites in RBD+NTD) > AUC of
  every baseline B1 variability-only, B2 exposure-only, B3 Hie-style 3-mer grammaticality +
  BLOSUM semantic change, B4 EVEscape-style independent-site frequency x exposure;
  10,000-draw bootstrap 95% CI on AUC difference vs best baseline excludes 0.
- Coordinate mapping SARS1->SARS2 spike via pairwise alignment of the two reference
  sequences (NC_045512 used ONLY as scoring coordinate frame), applied at scoring time.

## Case A predictions (SARS-CoV-1 spike numbering; mapped to SARS-CoV-2 only at scoring)
Top-30 RBD (318-510): 472, 489, 490, 432, 488, 486, 322, 476, 470, 445, 427, 480, 479, 487, 484, 446, 471, 485, 437, 442, 439, 428, 482, 443, 431, 468, 474, 392, 375, 391
Top-20 NTD (14-305): 206, 149, 198, 150, 294, 80, 70, 63, 37, 259, 154, 108, 33, 165, 215, 237, 166, 122, 211, 163

## Case B predictions (H3 mature HA numbering, HA1 1-329)
Top-30 HA1: 189, 144, 159, 137, 124, 158, 193, 156, 62, 133, 276, 225, 222, 226, 142, 145, 192, 155, 121, 278, 50, 186, 160, 135, 140, 143, 83, 126, 196, 188

## Predicted-site scores
RBD: 472:2.22; 489:1.87; 490:1.86; 432:1.84; 488:1.75; 486:1.72; 322:1.70; 476:1.69; 470:1.65; 445:1.65; 427:1.63; 480:1.62; 479:1.62; 487:1.60; 484:1.59; 446:1.59; 471:1.53; 485:1.51; 437:1.51; 442:1.49; 439:1.49; 428:1.48; 482:1.46; 443:1.46; 431:1.45; 468:1.39; 474:1.33; 392:1.30; 375:1.29; 391:1.28
NTD: 206:2.40; 149:2.26; 198:2.13; 150:2.04; 294:2.03; 80:2.00; 70:1.93; 63:1.93; 37:1.86; 259:1.86; 154:1.84; 108:1.83; 33:1.82; 165:1.82; 215:1.82; 237:1.81; 166:1.79; 122:1.79; 211:1.78; 163:1.77
HA1: 189:1.97; 144:1.93; 159:1.80; 137:1.73; 124:1.70; 158:1.66; 193:1.56; 156:1.54; 62:1.49; 133:1.47; 276:1.46; 225:1.44; 222:1.43; 226:1.42; 142:1.41; 145:1.40; 192:1.39; 155:1.34; 121:1.33; 278:1.30; 50:1.29; 186:1.27; 160:1.26; 135:1.25; 140:1.23; 143:1.21; 83:1.21; 126:1.21; 196:1.14; 188:1.14

SHA256(of body above) = d33a3f289bc0d882a32a82752713ab8488a1f16f9da3ed69d94fb73c297464ef
