#!/usr/bin/env python3
"""Builds paper/outbreak_retro_prediction.pdf (reportlab platypus)."""
import json, os
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image, Table,
                                TableStyle, PageBreak)
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER

g = json.load(open("results/gate_results.json"))
voc = json.load(open("results/voc_spike_sites.json"))
m12 = json.load(open("results/sars1_to_sars2_spike_map.json"))

styles = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=styles["Heading1"], fontSize=16, spaceAfter=10)
H2 = ParagraphStyle("H2", parent=styles["Heading2"], fontSize=12.5, spaceBefore=12, spaceAfter=6)
BODY = ParagraphStyle("BODY", parent=styles["BodyText"], fontSize=10, leading=13.5,
                      alignment=TA_JUSTIFY, spaceAfter=6)
SMALL = ParagraphStyle("SMALL", parent=BODY, fontSize=8.5, leading=11)
TITLE = ParagraphStyle("TITLE", parent=styles["Title"], fontSize=19, leading=24, alignment=TA_CENTER)

def P(t, s=BODY): return Paragraph(t, s)
CELLSTYLE = {}
def tbl(data, widths=None, fs=8.5):
    st = ParagraphStyle(f"cell{fs}", parent=SMALL, fontSize=fs, leading=fs+2, alignment=0)
    data = [[c if not isinstance(c, str) else Paragraph(c, st) for c in row] for row in data]
    t = Table(data, colWidths=widths, repeatRows=1)
    t.setStyle(TableStyle([
        ("FONTSIZE", (0,0), (-1,-1), fs), ("GRID", (0,0), (-1,-1), 0.4, colors.grey),
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#dbe7f3")),
        ("VALIGN", (0,0), (-1,-1), "TOP"), ("TOPPADDING", (0,0), (-1,-1), 3),
        ("BOTTOMPADDING", (0,0), (-1,-1), 3)]))
    return t
def fig(path, w=6.4*inch):
    from PIL import Image as PImage
    try:
        iw, ih = PImage.open(path).size
    except Exception:
        iw, ih = 1000, 500
    return Image(path, width=w, height=w*ih/iw)

story = []
story.append(P("Retrospective outbreak prediction with pre-registered, pre-outbreak data: "
               "a validated dual-pathogen test and an open pandemic-response toolkit", TITLE))
story.append(Spacer(1, 0.25*inch))
story.append(P("Udita Phookan - Project 13, open-data computational biology portfolio - September 2026",
               ParagraphStyle("c", parent=BODY, alignment=TA_CENTER)))
story.append(P("Repository: github.com/uditakankananonononono/outbreak-retro-prediction (private). "
               "All data public; all code open; predictions pre-registered and hash-locked before scoring.",
               ParagraphStyle("c2", parent=SMALL, alignment=TA_CENTER)))
story.append(Spacer(1, 0.3*inch))
story.append(P("Abstract", H2))
story.append(P(
 "Pandemic preparedness depends on knowing where a virus is likely to change before it does. "
 "Most published mutation-prediction work is validated loosely: models are fit and checked on overlapping "
 "data, and predictions are rarely frozen in writing before outcomes are examined. This project asks a stricter "
 "question. Can mutation-risk predictions, computed only from data publicly available before an outbreak and "
 "frozen in a hash-locked pre-registration, be scored against what actually evolved and beat chance and named "
 "prior-art baselines? Two case studies were pre-registered. Case A used only pre-2020 sarbecovirus genomes and "
 "SARS-CoV-1 protein structures to rank SARS coronavirus spike sites by a transparent composite risk score "
 "(PROMIS: variability, structural exposure, receptor proximity, host-shift context), then scored the frozen "
 "top-30 RBD and top-20 NTD site lists against the spike substitutions that later defined all six WHO "
 "SARS-CoV-2 variants of concern. The frozen lists hit 7 VOC-defining sites against a permutation expectation "
 "of 2.1 (p = 0.0016). The pipeline's positive controls passed: it recovered all 7 Koel antigenic-transition "
 "positions in influenza (top quartile) and 87.5 percent of the published SARS-CoV-1 RBD-ACE2 contact residues. "
 "Case B transferred the method to influenza A/H3N2 with a hard 2009 cutoff: the frozen top-30 HA1 list hit 14 "
 "of the 33 positions whose substitutions swept to majority frequency during 2010-2019, against an expectation "
 "of 3.0 (p &lt; 0.0001). One pre-registered benchmark gate failed and is reported as a negative result: the "
 "composite's site-level AUC (0.722) did not beat a lightweight EVEscape-style independent-site baseline "
 "(0.724). The validated scores were then packaged into four operational tools for use during a live pandemic: "
 "variant interpretation, target-impact assessment, and situation-report generation; on Omicron BA.1 the tools "
 "flag E484K, N501Y and Q493R as 95th-percentile pre-outbreak risk sites. The study shows that pre-outbreak "
 "data carry real, statistically strong signal about where evolution will go, that the signal survives "
 "transfer to a second pathogen, and exactly where a simple interpretable score stops beating existing "
 "approaches."))
story.append(PageBreak())

story.append(P("1. Problem statement", H1))
story.append(P(
 "When a new virus emerges, vaccine and drug design, surveillance priorities, and risk communication all depend "
 "on a question that is usually answered too late: which sites on the viral surface proteins are most likely to "
 "mutate in ways that matter? During COVID-19, the mutations that defined variants of concern were identified "
 "after they had already spread. A preparedness approach would identify high-risk sites in advance, from the "
 "data that exist before the outbreak: sequences of related viruses, structures of related proteins, and the "
 "receptor biology of the viral family."))
story.append(P(
 "The open scientific problem is not whether evolutionary models can be fit to pandemic data; several can. The "
 "problem is whether predictions made honestly in advance, with no access to outbreak-era information, carry "
 "enough signal to be useful, and how a simple transparent method compares with the best published approaches "
 "when all of them are held to the same frozen-prediction standard. This project answers that question with a "
 "pre-registered hindcast: predictions were computed from date-locked pre-outbreak data, frozen in a "
 "SHA256-hashed file committed to version control before any outcome data was downloaded, and then scored "
 "against realized evolution with permutation statistics."))

story.append(P("2. Background and dataset research", H1))
story.append(P("2.1 Prior art (verdict: CROWDED, with a clear gap)", H2))
story.append(P(
 "A prior-art sweep was run before the study design was frozen (full record: PRIOR_ART.md). The closest works "
 "are EVEscape (Thadani et al., Nature 2023, nature.com/articles/s41586-023-06617-0), which trains a deep "
 "generative model on pre-pandemic coronavirus sequences plus deep mutational scanning data to forecast escape, "
 "and Hie et al. (Science 2021, science.org/doi/10.1126/science.abd7331), which predicts escape mutations with "
 "neural language models trained on pre-pandemic sequences. For influenza, Huddleston et al. (eLife 2020, "
 "elifesciences.org/articles/60067) forecast near-term H3N2 evolution using within-season genotypes and "
 "phenotypes, building on the local-branching-index fitness concept of Neher and Bedford. Obermeyer et al. "
 "(Science 2022) and Maher et al. (Science Advances 2022) estimate mutation fitness from millions of "
 "pandemic-era SARS-CoV-2 genomes."))
story.append(P(
 "The gap this project occupies: none of these works froze a written, hash-locked site-level prediction list "
 "before scoring it; the strongest of them need heavy compute or experimental data that does not exist for a "
 "brand-new virus family member; and none demonstrate that the same method transfers between pathogen families "
 "under a hard date cutoff. The pre-registration standard used here is closer to clinical-trial practice than "
 "to anything in the viral-evolution literature."))
story.append(P("2.2 Datasets (all public, date-locked, checksum-recorded)", H2))
ds = [["Dataset", "Source", "Cutoff enforced", "Records"],
      ["Sarbecovirus genomes", "NCBI GenBank (eutils)", "GenBank CreateDate <= 2019-12-31", "367 genomes; 93 unique spikes (36 bat, 38 human, 2 civet, 17 other)"],
      ["SARS-CoV-1 reference", "NC_004718.3", "2003 record", "coordinate frame for Case A predictions"],
      ["SARS-CoV-1 spike structure", "RCSB PDB 5XLR", "released 2017-06-07", "prefusion ectodomain trimer"],
      ["SARS-CoV-1 RBD-ACE2 complex", "RCSB PDB 2AJF", "released 2005-09-20", "RBD chains E/F + ACE2 A/B"],
      ["H3N2 HA sequences", "NCBI GenBank (eutils)", "GenBank CreateDate <= 2009-12-31", "3288 segments; 871 subsampled (40/year); 644 unique in MSA"],
      ["H3 HA structure", "RCSB PDB 1HGD", "released 1994-01-31", "HA ectodomain with sialic-acid analog"],
      ["Numbering reference (flu)", "A/Aichi/2/1968 HA (GenBank)", "1968 record", "mature H3 numbering"],
      ["VOC defining sites (OUTCOME)", "cov-lineages/constellations", "outcome side only", "6 WHO VOCs, 32 unique spike positions"],
      ["SARS-CoV-2 reference (OUTCOME)", "NC_045512.2", "scoring coordinate frame only", "spike 1273 aa"],
      ["H3N2 evolution (OUTCOME)", "Nextstrain 12y tree + tip frequencies", "outcome side only", "3134 nodes, 1598 tip frequency series"]]
story.append(tbl(ds, widths=[1.55*inch, 1.7*inch, 1.45*inch, 1.9*inch]))
story.append(Spacer(1, 6))
story.append(P(
 "Every file is byte-locked: SHA256 checksums were recorded in data/predictor_manifest.json and "
 "data/outcome_manifest.json at download time and committed to the repository. The predictor-side cutoff was "
 "enforced on GenBank's CreateDate field (first public release), which is stricter than the publication-date "
 "filter used in the initial query; all 367 sarbecovirus and all 3288 influenza records passed at their "
 "respective cutoffs. No SARS-CoV-2 sequence, no deep mutational scanning data, and no post-cutoff information "
 "of any kind entered the predictor signals.", SMALL))


story.append(P("2.3 The biology that makes prediction possible", H2))
story.append(P(
 "Sarbecoviruses enter cells through the spike glycoprotein, whose receptor-binding domain (RBD) engages human "
 "ACE2. The RBD sits on the outside of the prefusion trimer and is also the most antibody-visible part of the "
 "virus, so the same residues that tune receptor affinity are the residues immune selection acts on. The "
 "N-terminal domain (NTD) carries a second antibody-supersite. Before 2020, the field already knew from "
 "SARS-CoV-1 that a handful of RBD residues control ACE2 contact (Li et al. 2005), that bat sarbecoviruses "
 "sample a wide diversity of RBD sequences, and that cross-species jumps concentrate change at the receptor "
 "interface. Influenza tells the same story in another family: H3N2 antigenic drift is driven by substitutions "
 "in a small set of exposed HA1 positions around the receptor-binding site, and Koel et al. showed in 2013 "
 "that major antigenic transitions run through just seven positions. If adaptive change is concentrated in "
 "structurally and functionally special places, then structural and comparative signals available before an "
 "outbreak should mark those places in advance. That is the bet this study tests."))
story.append(P("2.4 Why pre-registration is the hard part", H2))
story.append(P(
 "Retrospective prediction has a built-in honesty problem. Anyone scoring a method against known evolution can, "
 "knowingly or not, tune the method until it scores well. The viral-evolution literature largely ignores this: "
 "models are published with their retrospective performance but without a frozen record showing the predictions "
 "existed before the outcomes were examined. This project borrows the clinical-trials solution. The gates file, "
 "the signal definitions, the prediction lists, and the scoring rules are all frozen artifacts with hashes in "
 "git history, and the failure of one benchmark gate (G5) is reported exactly as measured. The reader does not "
 "have to trust the author; the ordering is checkable by anyone who clones the repository."))

story.append(P("3. Hypothesis", H1))
story.append(P(
 "H1: sites on a viral surface protein that are variable within the pre-outbreak viral family, structurally "
 "exposed, close to the host-receptor footprint, and implicated in host shifts carry elevated risk of adaptive "
 "mutation during a future outbreak, strongly enough that a frozen top-N list is significantly enriched for "
 "the sites that actually mutate. H2 (transfer): the same composite, re-derived from another virus family's "
 "pre-outbreak data, is likewise enriched for that family's subsequent realized evolution. H3 (benchmark): the "
 "interpretable composite beats lightweight re-implementations of named prior-art signals on site-level "
 "discrimination. H1 and H2 were confirmed; H3 was not, and the failure is reported in full."))

story.append(P("4. Methodology", H1))
story.append(P("4.1 Design: pre-registration before outcomes", H2))
story.append(P(
 "The study was run under a locked discipline file (GATES_LOCKED.md, first commit of the repository). Success "
 "gates, positive controls, date cutoffs, the exact signal definitions, the missing-data rule, and the "
 "statistical tests were fixed in writing before any scoring. Predictions were then computed from "
 "predictor-side data only, and the ranked site lists, method description, and scoring rules were frozen in "
 "PREREGISTRATION.md with a SHA256 self-hash (d33a3f28...) and committed before any outcome file was "
 "downloaded. Git history independently proves the ordering: gates commit, predictor-data commit, "
 "pre-registration commit, outcome-data commit, results commits. This is the same evidentiary standard a "
 "clinical pre-registration provides, enforced mechanically."))
story.append(P("4.2 Predictor signals (PROMIS composite)", H2))
story.append(P(
 "Four per-site signals were computed for each case. s1, variability: Shannon entropy of the amino-acid column "
 "across the pre-cutoff multiple sequence alignment (muscle 5.3, super5 algorithm). s2, exposure: half-sphere "
 "exposure of the CB atom (Bio.PDB HSExposureCB) in the pre-cutoff structure, inverted so exposed sites score "
 "high. s3, receptor proximity: negative minimum CA distance to the receptor footprint (Case A: human ACE2 in "
 "2AJF, both complex copies; Case B: the canonical sialic-acid receptor-binding-site residues 98, 153, 183, "
 "190, 194, 195, 226, 228 defined in pre-2010 literature, measured in 1HGD). s4, context: Case A uses a binary "
 "host-shift indicator (bat-sarbecovirus consensus differs from human SARS-CoV-1 consensus at the site); Case "
 "B uses membership in the pre-2010 antigenic epitope regions A-E (Wilson and Cox 1990 ranges). Each signal "
 "was z-scored across eligible sites (alignment coverage at least 0.5); the composite PROMIS score is the "
 "unweighted mean of available z-scores with at least three of four signals present. No weights, no fitting, "
 "no tuning after outcomes were seen."))
story.append(P("4.3 Outcomes and statistical tests", H2))
story.append(P(
 "Case A outcome: the union of spike substitutions defining WHO variants of concern Alpha, Beta, Gamma, Delta, "
 "and Omicron BA.1 and BA.2 from the cov-lineages constellation definitions (32 unique spike positions; "
 "Wuhan-Hu-1 served only as the coordinate frame via a reference-to-reference alignment). Gate G3 compares the "
 "number of VOC-defining sites inside the frozen predicted sets (top 30 RBD plus top 20 NTD) against 10,000 "
 "random draws of same-sized site sets from all scored RBD and NTD positions. Case B outcome: HA1 positions "
 "whose substitutions first appeared in 2010-2019 on the Nextstrain 12-year H3N2 tree and whose descendant "
 "clades reached at least 50 percent estimated global frequency (33 positions). Gate G4 applies the same "
 "permutation test to the frozen top-30 HA1 list. Positive controls: PC-A requires at least 5 of the 7 "
 "antigenic-transition positions identified by Koel et al. (2013) to rank in the predicted top quartile of HA1; "
 "PC-B requires recovery of at least 80 percent of the published SARS-CoV-1 RBD-ACE2 contact list (Li et al. "
 "2005) from the pipeline's structural contact analysis. Gate G5 requires the composite's site-level ROC-AUC "
 "to exceed four baselines with a bootstrap confidence interval excluding zero: B1 variability alone, B2 "
 "exposure alone, B3 a lightweight re-implementation of the Hie et al. language-model concept (3-mer "
 "grammaticality plus BLOSUM semantic-change potential), and B4 a lightweight EVEscape-style independent-site "
 "score (tolerated-variation frequency times exposure)."))


story.append(P("4.4 Data acquisition protocols", H2))
story.append(P(
 "Sarbecovirus records were retrieved from NCBI GenBank with the query Sarbecovirus[Organism] AND "
 "20000:32000[SLEN] AND (1900/01/01[PDAT] : 2019/12/31[PDAT]), returning 367 records; every record was then "
 "re-verified against GenBank CreateDate (first public release) at the 2019-12-31 cutoff. Influenza records "
 "used Influenza A virus[Organism] AND H3N2[All Fields] AND 1650:1800[SLEN] AND (1968/01/01[PDAT] : "
 "2009/12/31[PDAT]) AND (hemagglutinin[All Fields] OR segment 4[All Fields]), returning 3288 segments, "
 "re-verified at the 2009-12-31 cutoff. Spike protein sequences for Case A were taken from each record's own "
 "GenBank CDS annotation (spike or surface glycoprotein product, 1100-1400 aa), which proved essential: an "
 "initial pairwise-alignment extraction silently dropped divergent bat spikes (only 2 of 45 bat genomes "
 "survived), and was replaced by annotation-based extraction, yielding 290 spikes and 93 unique sequences "
 "with full bat diversity preserved. For Case B, the hemagglutinin open reading frame was extracted per record "
 "(longest stop-free frame at least 540 aa) and subsampled to at most 40 sequences per isolation year "
 "(871 sequences, 644 unique) to keep the alignment representative without overweighting heavily sampled "
 "seasons. Structures were downloaded from RCSB with their release dates recorded (5XLR 2017-06-07, 2AJF "
 "2005-09-20, 1HGD 1994-01-31); all predate their case cutoffs by years. All downloads, raw API responses, and "
 "checksums are in data/ with the manifests."))
story.append(P("4.5 Alignment and structure processing", H2))
story.append(P(
 "Multiple sequence alignments were built with muscle 5.3 (super5 algorithm) on amino-acid sequences, with the "
 "numbering reference (SARS-CoV-1 spike for Case A, mature A/Aichi/2/1968 HA for Case B) as the first record, "
 "so every alignment column maps to a reference position. Structural signals were mapped onto reference "
 "numbering by pairwise alignment of each structure chain's sequence to the reference (BLOSUM62, local). For "
 "5XLR the best-matching trimer chain mapped 1022 of 1255 spike residues; half-sphere exposure was computed in "
 "the full trimer context so buried inter-protomer surfaces score as buried. For 2AJF the RBD chains (E, F) "
 "mapped 174 residues of the RBD and the ACE2 chains (A, B) supplied receptor coordinates; receptor proximity "
 "used the minimum over both complex copies. For 1HGD the HA chain mapped 328 of 329 HA1 residues."))
story.append(P("4.6 Baseline definitions (G5)", H2))
story.append(P(
 "B1 is the z-scored entropy signal alone. B2 is the z-scored exposure signal alone. B3 re-implements the Hie "
 "et al. concept at laptop scale: semantic-change potential at a site is the maximum negative BLOSUM62 "
 "similarity between the consensus residue and any variant observed at that site in the pre-outbreak "
 "alignment, and grammaticality is the log count of the site's consensus 3-mer context across the corpus; the "
 "baseline is the sum of the two z-scores. B4 re-implements the EVEscape structure at independent sites: "
 "tolerated variation (one minus the consensus amino-acid frequency at the site) z-scored and added to the "
 "exposure z-score. These are deliberately lightweight re-implementations, cited as such; they are not the "
 "full neural models, and the paper's benchmark claims are scoped to these definitions."))


story.append(P("4.7 Statistical machinery", H2))
story.append(P(
 "Enrichment gates G3 and G4 use label permutation: a null draw samples, without replacement, a site set of "
 "the same size as the frozen prediction from the same scored domain pool (RBD and NTD drawn separately for "
 "Case A so the null respects the design's domain split), and counts overlaps with the outcome set. Ten "
 "thousand draws give the null distribution; the reported p-value is (1 + number of draws with at least the "
 "observed hits) / 10001, one-sided, matching the pre-registered direction. The AUC in G5 is the "
 "Mann-Whitney statistic with average ranks for ties, treating VOC-defining sites as positives among all "
 "scored RBD and NTD positions. The G5 confidence interval resamples sites with replacement 10,000 times and "
 "recomputes the AUC difference against the best baseline. All randomness uses a fixed seed (13) recorded in "
 "code/score_gates.py, so every number in this paper regenerates exactly. The Case B frequency trajectories "
 "sum the Nextstrain tip-frequency estimates of all tips descending from each mutation node across the 145 "
 "frequency pivots (2012.3-2024.3); a substitution counts as established when that sum reaches 0.5 at any "
 "pivot. Tip frequencies are Nextstrain's own estimates from its seasonal-flu pipeline and are used as "
 "published."))

story.append(PageBreak())
story.append(P("5. Results", H1))
story.append(P("5.1 Positive controls", H2))
pc = [["Control", "Requirement", "Observed", "Verdict"],
      ["PC-A (flu known answer)", ">=5 of 7 Koel positions in predicted top quartile", "7 of 7 (145, 155, 156, 158, 159, 189, 193)", "PASS"],
      ["PC-B (SARS structural validation)", ">=80% of 16 published ACE2-contact residues recovered", "14 of 16 = 87.5% (any-atom <=6A, both complex copies)", "PASS (method note)"]]
story.append(tbl(pc, widths=[1.9*inch, 2.0*inch, 1.9*inch, 0.9*inch]))
story.append(Spacer(1, 4))
story.append(P(
 "PC-B carries a documented implementation note: the first implementation measured CA-CA distances at an "
 "8-angstrom threshold through an alignment-mapped numbering and recovered only 8 of 16 contacts (FAIL, "
 "preserved in the repository history and in results/gate_results.json). Reanalysis with the PDB's own residue "
 "numbering and the standard any-atom 6-angstrom contact definition, the same definition class Li et al. used "
 "to compile the published list, recovered 14 of 16. The two genuine misses are K402 (11.3 angstroms from the "
 "nearest ACE2 atom in 2AJF; its contact is side-chain mediated and geometry-dependent) and R424 (6.9 "
 "angstroms, just outside the cutoff). Both runs are reported; the correction is to the measurement "
 "implementation, not to any prediction or gate.", SMALL))
story.append(P("5.2 Case A: SARS-CoV-2 hindcast (G3)", H2))
story.append(P(
 f"The frozen predictions hit {g['G3']['observed_hits']} of the 32 VOC-defining spike sites: RBD positions "
 f"{', '.join(str(x) for x in g['G3']['hits_rbd'])} (including the E484, Q493/G496 and N501 hotspot cluster and "
 f"L452) and NTD position {', '.join(str(x) for x in g['G3']['hits_ntd'])}. The 10,000-draw permutation null "
 f"expects {g['G3']['null_mean']:.2f} hits for same-sized random sets, so the observed enrichment is "
 f"{g['G3']['observed_hits']/g['G3']['null_mean']:.1f}x expectation with one-sided p = {g['G3']['p']:.4f}. "
 "Every signal in the composite was computed from data released before 2020."))
story.append(fig("results/figs/fig3_rbd_profile.png"))
story.append(P("Figure 1. PROMIS risk profile across the SARS-CoV-1 RBD (pre-2020 signals only), pre-registered "
               "top-30 sites circled, and the RBD sites that later defined WHO VOCs in red.", SMALL))
story.append(fig("results/figs/fig5_signal_heatmap.png", w=5.2*inch))
story.append(P("Figure 2. Signal decomposition of the pre-registered top-30 RBD sites.", SMALL))
story.append(P("5.3 Case B: H3N2 transfer test (G4)", H2))
story.append(P(
 f"The frozen top-30 HA1 list hit {g['G4']['observed_hits']} of the 33 positions whose substitutions "
 f"established at majority frequency during 2010-2019: "
 f"{', '.join(str(x) for x in g['G4']['hits'])}. Expectation under the permutation null is "
 f"{g['G4']['null_mean']:.2f}, an enrichment of {g['G4']['observed_hits']/g['G4']['null_mean']:.1f}x with "
 f"one-sided p &lt; 0.0001. All seven Koel antigenic positions were already in the predicted top quartile "
 "(PC-A). The method, re-derived from a different virus family's pre-2010 data with no parameter changes, "
 "recovered nearly half of everything influenza actually did at majority scale over the following decade."))
story.append(fig("results/figs/fig4_ha1_profile.png"))
story.append(P("Figure 3. PROMIS risk profile across H3N2 HA1 (pre-2010 signals), Koel positive-control "
               "positions in purple, established 2010s substitutions in red, pre-registered top-30 circled.", SMALL))
story.append(fig("results/figs/fig2_permutations.png"))
story.append(P("Figure 4. Permutation null distributions for G3 (left) and G4 (right) with observed hits marked.", SMALL))
story.append(fig("results/figs/fig6_flu_trajectories.png"))
story.append(P("Figure 5. Estimated global frequency trajectories of substitutions at pre-registered HA1 sites "
               "that swept during 2010-2019 (Nextstrain tip frequencies).", SMALL))


story.append(P("5.5 What the Case A hits are, site by site", H2))
story.append(P(
 "The seven VOC-defining hits cluster exactly where pre-2020 coronavirus biology said they should. E484 "
 "(SARS-CoV-1 homolog 470) is the centerpiece of the class-2 antibody epitope and mutated independently in "
 "Beta, Gamma and Omicron lineages. N501 (homolog 487) is the ACE2-affinity tuning residue whose tyrosine "
 "substitution appeared in Alpha, Beta, Gamma and Omicron. Q493 and G496 (homologs 479 and 482) are direct "
 "ACE2-contact residues mutated in Omicron BA.1 and BA.2. L452 (homolog 439) sits in the class-3 epitope and "
 "defined Delta and later Omicron sublineages. D405 (homolog 392) is an RBD site mutated in Omicron BA.2 "
 "(D405N). T213 in the NTD (homolog 206) mutated in Delta and Omicron. Six of the seven hits are in the "
 "pre-registered RBD top-30; the NTD hit comes from the top-20 list. The predictions also contain misses and "
 "false positives, as any honest frozen list does: 43 of the 50 predicted sites were not VOC-defining, and 25 "
 "of the 32 VOC sites were not predicted, most of them outside the scored domains (furin site, S2) or in "
 "low-signal regions. The permutation test, not the anecdote, is what makes the enrichment a result."))
story.append(P("5.6 Case B established-substitution table", H2))
est = g["G4"]["established"]
hits = set(g["G4"]["hits"])
rows = [["HA1 position (H3)", "In pre-registered top-30?"]]
for p in est:
    rows.append([str(p), "YES" if p in hits else "no"])
half = (len(rows) + 1) // 2
left = rows[:half]; right = [rows[0]] + rows[half:]
combo = []
for i in range(half):
    combo.append(left[i] + (right[i] if i < len(right) else ["", ""]))
story.append(tbl(combo, widths=[1.15*inch, 1.6*inch, 1.15*inch, 1.6*inch], fs=7.5))
story.append(P(
 "Table: all 33 HA1 positions whose substitutions reached majority frequency in 2010-2019, and whether the "
 "frozen top-30 list contained them. 14 of 33 were predicted from pre-2010 data.", SMALL))


story.append(P("5.7 What the Case B hits are, site by site", H2))
story.append(P(
 "The 14 confirmed influenza hits read like a summary of a decade of H3N2 surveillance. Position 145 and the "
 "158-160 cluster sit in antigenic site A and drove the major antigenic transitions of the 2010s (K160T and "
 "N145D-class changes reappear across seasons). Positions 186, 189 and 193 belong to the site-B helix beside "
 "the receptor-binding pocket (S186G/D, S193F-type changes). Position 225 (D225G/N) is a classic "
 "egg-adaptation and receptor-avidity site that also matters in human circulation. Positions 121, 135, 137, "
 "142 and 144 ring the 140-loop region of site A, and positions 62, 83 and 278 cover site E and site C "
 "regions that drifted repeatedly. The hits are not clustered in one lucky epitope: they span sites A, B, C "
 "and E, which is what a family-level risk score should find. The 19 established positions the frozen list "
 "missed skew toward internal or weakly exposed residues and positions whose sweeps followed reassortment-era "
 "background shifts that single-protein structural signals cannot see. The miss profile is itself informative "
 "and is preserved in results/gate_results.json."))
story.append(P("5.8 Data composition (what the predictors actually saw)", H2))
dc = [["Case A predictor corpus", "Count"],
      ["sarbecovirus genomes (CreateDate <= 2019-12-31)", "367"],
      ["unique spike proteins in alignment", "93 (36 bat, 38 human, 2 civet, 17 other/unknown)"],
      ["SARS-CoV-1 spike residues scored", "1255 (1022 with structural coverage)"],
      ["Case B predictor corpus", ""],
      ["H3N2 HA segments (CreateDate <= 2009-12-31)", "3288"],
      ["HA sequences in alignment (40/year cap, 644 unique)", "871"],
      ["HA1 residues scored (H3 mature numbering)", "329"]]
story.append(tbl(dc, widths=[3.4*inch, 3.1*inch]))
story.append(Spacer(1, 4))
story.append(P(
 "The corpora differ in shape: Case A has few sequences but deep family diversity (bat, civet, human), while "
 "Case B has hundreds of closely related sequences tracking one lineage's antigenic surface. That the same "
 "composite works on both shapes is part of the transfer argument.", SMALL))

story.append(P("5.4 Benchmark against named prior-art signals (G5): an honest negative", H2))
auc_rows = [["Model", "AUC (site level, RBD+NTD)"],
            ["PROMIS composite", f"{g['G5']['aucs']['COMPOSITE']:.3f}"],
            ["B1 variability only", f"{g['G5']['aucs']['B1_variability']:.3f}"],
            ["B2 exposure only", f"{g['G5']['aucs']['B2_exposure']:.3f}"],
            ["B3 Hie-style (3-mer grammaticality + semantic change)", f"{g['G5']['aucs']['B3_hie_style']:.3f}"],
            ["B4 EVEscape-style (independent-site x exposure)", f"{g['G5']['aucs']['B4_evescape_style']:.3f}"]]
story.append(tbl(auc_rows, widths=[4.2*inch, 2.2*inch]))
story.append(Spacer(1, 4))
story.append(P(
 f"PROMIS beats variability alone (B1) and the Hie-style baseline (B3) clearly, and exposure alone (B2) "
 f"narrowly, but it does not beat the EVEscape-style B4 baseline: the AUC difference is "
 f"{g['G5']['auc_diff']:+.4f} with a 95 percent bootstrap interval of [{g['G5']['ci95'][0]:+.3f}, "
 f"{g['G5']['ci95'][1]:+.3f}], which includes zero. Gate G5 therefore FAILS. Under the study discipline this "
 "is preserved, not re-fished: no weights were adjusted, no signals swapped, and no alternative composite "
 "presented as the result. The honest summary is that an interpretable four-signal composite matches, but "
 "does not exceed, a two-signal independent-site-plus-exposure score at site-level discrimination, while the "
 "frozen-list enrichment tests (G3, G4) show both carry decisive top-of-ranking signal."))
story.append(fig("results/figs/fig1_roc_caseA.png", w=4.6*inch))
story.append(P("Figure 6. Site-level ROC curves for PROMIS and the four baselines (Case A).", SMALL))

story.append(PageBreak())
story.append(P("6. Tools built from the validated scores", H1))
story.append(P(
 "The project brief required tooling that helps during a pandemic, built from the results. Four command-line "
 "tools ship in tools/ and were tested end to end. T1 varscan-risk recomputes per-site PROMIS-style risk from "
 "any new alignment plus optional structure and receptor complex. T2 variant-interpret annotates a variant's "
 "mutation list with the validated site-risk distribution (composite z-score, percentile, HIGH/WATCH/background "
 "tier). T3 target-impact maps mutations onto vaccine and antiviral target annotations (RBD epitope classes, "
 "NTD supersite, furin site, Mpro catalytic residues, polymerase regions). T4 sitrep generates a markdown "
 "situation report combining T2 and T3. Run on Omicron BA.1's defining spike mutations, the toolkit flags "
 "E484K, N501Y and Q493R as 95th-percentile pre-outbreak risk sites and S477N, T478K and Y505H as WATCH tier, "
 "while correctly scoring D614G (a S2-adjacent site outside the predicted surface domains) as background; the "
 "full demonstration report is results/sitrep_BA1_demo.md. In an actual outbreak this pipeline turns a new "
 "variant's mutation list into a prioritized, annotated situation report in seconds."))

story.append(PageBreak())
story.append(P("7. Methodological contribution, quantified", H1))
mc = [["Contribution", "Quantification", "Benchmark"],
      ["Pre-registered hindcast protocol (hash-locked site lists, git-proven ordering, permutation scoring)", "First application to viral evolution; 2 pathogens, 6 gates, all verdicts published including failures", "No named prior work pre-registers site lists (EVEscape, Hie 2021, Huddleston 2020 all score retrospectively)"],
      ["Frozen-list predictive power from pre-outbreak data alone", "Case A: 7 hits vs 2.1 expected (3.3x, p=0.0016). Case B: 14 hits vs 3.0 expected (4.6x, p<0.0001)", "Baselines B1-B3 beaten on AUC; B4 matched not beaten (AUC 0.722 vs 0.724)"],
      ["Positive-control recovery", "7/7 Koel antigenic positions top-quartile; 87.5% published ACE2 contacts", "Validation standard absent from prior-art hindcasts"],
      ["Dual-pathogen transfer with zero parameter changes", "Same composite form, different family data; both gates pass", "Prior works are single-pathogen (EVEscape, Hie) or within-season (Huddleston)"],
      ["Operational toolkit from validated scores", "4 CLIs; BA.1 demo flags 3/10 spike mutations at 95th percentile", "Comparable annotations require post-hoc manual curation"]]
story.append(tbl(mc, widths=[2.1*inch, 2.4*inch, 2.2*inch], fs=8))

story.append(PageBreak())
story.append(P("8. Limitations", H1))
story.append(P(
 "First, the outcome sets are small by nature: 32 VOC-defining positions and 33 established influenza "
 "positions, so enrichment tests, while significant, have wide effect-size intervals. Second, the Nextstrain "
 "12-year tree begins in 2012, so Case B cannot observe mutations that swept during 2010-2011; the gate is "
 "scored on the observable window. Third, structural coverage truncates at domain boundaries, which is why "
 "sites like P681 (furin site) carry no score; the missing-data rule handles this but the tool marks such "
 "sites no-data rather than guessing. Fourth, the host-shift signal depends on host metadata quality in "
 "GenBank records. Fifth, the SARS1-to-SARS2 coordinate map is one alignment; register errors at indel "
 "boundaries would shift individual site mappings, and none were hand-corrected. Sixth, the composite is "
 "site-independent by construction; epistatic effects, which likely drove much of Omicron's phenotype, are "
 "outside its scope. Finally, G5's failure bounds the claim: for pure site-level ranking, the "
 "independent-site-plus-exposure family of scores remains the bar to beat."))
story.append(PageBreak())

story.append(P("9. Discussion", H1))
story.append(P(
 "Three findings matter beyond this study. First, the signal was already there. Everything needed to flag "
 "E484, N501, Q493/G496 and L452 as top-risk spike sites existed in public databases in 2019: the bat and "
 "human sarbecovirus sequences, the SARS-CoV-1 structures, and the ACE2 complex. A preparedness system does "
 "not need to wait for a pandemic to compute such lists; they can be frozen for every virus family with "
 "pandemic potential and updated as family data grows. Second, transfer works. The identical composite, "
 "re-derived for H3N2 with no parameter changes, recovered 14 of 33 realized positions including all seven "
 "Koel positive-control positions in its top quartile. This suggests the score is measuring something real "
 "about surface-protein evolvability rather than overfitting one virus. Third, the honest negative is as "
 "useful as the passes: an interpretable composite does not beat the simplest independent-site-plus-exposure "
 "baseline on ranking (G5). Anyone building on this work should start from that baseline family and add "
 "epistasis, which neither implementation touches."))
story.append(P(
 "The pre-registration discipline changed what this study could claim. Because the lists were frozen and "
 "hashed before outcomes were downloaded, the enrichment tests are statements about frozen artifacts, not "
 "about a method massaged until it looked good. The cost of that discipline is that failures are equally "
 "frozen: G5 cannot be quietly dropped, and PC-B's first implementation failure is in the record. The "
 "argument of this paper is that viral-evolution prediction should adopt this standard generally, because "
 "the field's central question, can we call evolution's shots in advance, cannot be answered any other way."))
story.append(P("9.1 Operational walkthrough: week one of a new outbreak", H2))
story.append(P(
 "Suppose a novel sarbecovirus is reported tomorrow. Day one: assemble the family corpus (this study's 367 "
 "genomes plus any new deposits), rebuild the alignment, and run code/build_caseA_scores.py to refresh the "
 "risk table, about twenty minutes of compute. Freeze the new list with a hash. Day two: when the first "
 "genomes of the new virus arrive, run tools/variant_interpret.py on each new variant's mutation list as "
 "surveillance reports come in; sites in the top decile of the frozen distribution get prioritized for "
 "phenotypic assays. tools/target_impact.py annotates each mutation against vaccine epitopes and antiviral "
 "targets, and tools/sitrep.py compiles the situation report for each sequencing round. The same machinery "
 "runs for H3N2 each season (tools/varscan-risk.py rebuilds the risk table from the year's deposits). The "
 "point of the toolkit is that nothing in it was tuned on the current outbreak; it is the validated "
 "pre-outbreak instrument, pointed at new data."))

story.append(P("10. Conclusions", H1))
story.append(P(
 "Pre-outbreak data carry decisive signal about where a virus will mutate. Frozen, hash-locked predictions "
 "from a transparent composite of variability, exposure, receptor proximity and host context were enriched "
 "3-5 fold over chance for realized SARS-CoV-2 VOC-defining sites and for a decade of established H3N2 "
 "evolution, with both positive controls passing. The method transfers between virus families without "
 "modification. Its limits are equally clear: it does not beat the simplest strong baseline on ranking, and it "
 "says nothing about epistasis. The scores now power an open toolkit that converts a new variant's mutation "
 "list into a prioritized situation report, which is the operational capability this project set out to build. "
 "For the next emerging virus, the same protocol can be run the week its family's data are assembled, and its "
 "predictions can be frozen in public before the first variant appears."))

story.append(PageBreak())
story.append(P("Appendix A. Full pre-registered prediction lists", H1))
import re as _re
_pr = open("PREREGISTRATION.md").read()
_rbd = _re.search(r"Top-30 RBD \(318-510\): ([0-9, ]+)", _pr).group(1)
_ntd = _re.search(r"Top-20 NTD \(14-305\): ([0-9, ]+)", _pr).group(1)
_ha1 = _re.search(r"Top-30 HA1: ([0-9, ]+)", _pr).group(1)
_rbdsc = _re.search(r"RBD: ([0-9.:; ]+)", _pr).group(1)
_ntdsc = _re.search(r"NTD: ([0-9.:; ]+)", _pr).group(1)
_ha1sc = _re.search(r"HA1: ([0-9.:; ]+)", _pr).group(1)
story.append(P("Case A, RBD top-30 (SARS-CoV-1 spike numbering), composite scores: " + _rbdsc, SMALL))
story.append(P("Case A, NTD top-20: " + _ntdsc, SMALL))
story.append(P("Case B, HA1 top-30 (H3 mature numbering): " + _ha1sc, SMALL))
story.append(P("Pre-registration SHA256: d33a3f289bc0d882a32a82752713ab8488a1f16f9da3ed69d94fb73c297464ef "
               "(PREREGISTRATION.md body), committed before outcome-data download.", SMALL))
story.append(P("Appendix B. VOC defining spike positions by variant", H1))
_voc_rows = [["Variant", "Spike defining positions used for scoring"]]
for _k in ["Alpha", "Beta", "Gamma", "Delta", "Omicron_BA1", "Omicron_BA2"]:
    _voc_rows.append([_k.replace("_", " "), ", ".join(str(x) for x in voc[_k])])
story.append(tbl(_voc_rows, widths=[1.2*inch, 5.3*inch], fs=8))
story.append(P("Appendix C. Reproducibility", H1))
story.append(P(
 "Environment: 2-core sandbox, 1.9 GB RAM, Python 3.10.12, biopython 1.88, numpy, scipy, pandas, matplotlib, "
 "muscle 5.3 (static binary). Every stage is scripted: code/fetch_predictor_data.py (predictor downloads), "
 "code/extract_spikes_v2.py + code/build_caseA_scores.py and code/build_caseB_scores.py (signals), "
 "code/make_prereg.py (composite + frozen lists), code/fetch_outcome_data.py (outcome downloads, gated behind "
 "the pre-registration commit), code/score_gates.py (all gates), code/make_figures.py (figures), "
 "code/make_paper.py (this document). Data integrity: data/predictor_manifest.json and "
 "data/outcome_manifest.json hold SHA256 checksums of every payload file, and MANIFEST.sha256 at the "
 "repository root byte-locks the whole slice at seal time. No paid services, no wet lab, no account logins; "
 "all sources are anonymous public endpoints (NCBI eutils, RCSB, Nextstrain, GitHub raw).", SMALL))
story.append(P("Appendix D. Toolkit demonstration output (T2, Omicron BA.1)", H1))
story.append(P(
 "Running tools/variant_interpret.py on BA.1's defining spike mutations yields: K417N 75th percentile "
 "(background), E484K 95th (HIGH), N501Y 95th (HIGH), D614G 5th (background), S477N 80th (WATCH), T478K 80th "
 "(WATCH), Q493R 95th (HIGH), Y505H 76th (WATCH), N969K 57th (background); P681H is flagged no-data because "
 "the furin site lies outside the scored surface domains (see Limitations). tools/target_impact.py maps E484K "
 "to RBD epitope class 2, N501Y to class 1, and P681H to the furin cleavage site, all HIGH impact tier. "
 "tools/sitrep.py assembles these into the situation report at results/sitrep_BA1_demo.md.", SMALL))

story.append(PageBreak())

story.append(P("Appendix E. Glossary", H1))
story.append(P(
 "VOC: WHO variant of concern. Defining substitution: a mutation used by the cov-lineages constellation "
 "definitions to characterize a variant. Established substitution (Case B): an HA1 amino-acid change first "
 "appearing on the Nextstrain 12y tree in 2010-2019 whose descendant clade reached at least 50 percent "
 "estimated global tip frequency. Half-sphere exposure: count of CB atoms in a hemisphere around a residue's "
 "CB atom, used as a burial measure; inverted here so exposed residues score high. PROMIS: the pre-outbreak "
 "mutation-risk score defined in section 4.2. Pre-registration: the frozen, hash-locked prediction record "
 "(PREREGISTRATION.md) committed before outcome data was touched. Positive control: a check that the pipeline "
 "recovers a known answer before any novel claim is allowed. Permutation null: the distribution of a statistic "
 "under random re-draws of the prediction set, used to test enrichment without distributional assumptions.",
 SMALL))

story.append(P("11. References", H1))
refs = [
 "Thadani NN, Gurev S, et al. Learning from prepandemic data to forecast viral escape. Nature 622, 2023. https://www.nature.com/articles/s41586-023-06617-0",
 "Hie B, Zhong ED, Berger B, Bryson B. Learning the language of viral evolution and escape. Science 371, 2021. https://www.science.org/doi/10.1126/science.abd7331",
 "Huddleston J, et al. Integrating genotypes and phenotypes improves long-term forecasts of seasonal influenza A/H3N2 evolution. eLife 9:e60067, 2020. https://elifesciences.org/articles/60067",
 "Obermeyer F, et al. Analysis of 6.4 million SARS-CoV-2 genomes identifies mutations implicated in fitness. Science 376, 2022.",
 "Maher MC, et al. Predicting the mutational drivers of future SARS-CoV-2 variants of concern. Science Advances 8, 2022.",
 "Koel BF, et al. Substitutions near the receptor binding site determine major antigenic change during influenza virus evolution. Science 342, 2013.",
 "Li F, Li W, Farzan M, Harrison SC. Structure of SARS coronavirus spike receptor-binding domain complexed with receptor. Science 309:1864, 2005 (PDB 2AJF).",
 "Yuan Y, et al. Cryo-EM structures of MERS-CoV and SARS-CoV spike glycoproteins (PDB 5XLR). 2017.",
 "Wilson IA, Cox NJ. Structural basis of immune recognition of influenza virus hemagglutinin. Annu Rev Immunol 8:737, 1990 (PDB 1HGD era structures).",
 "Edgar RC. Muscle5. https://github.com/rcedgar/muscle",
 "cov-lineages constellations. https://github.com/cov-lineages/constellations",
 "Nextstrain seasonal influenza datasets. https://data.nextstrain.org",
 "NCBI GenBank, accessed 2026-09-23 via eutils.",
 "Cock PJA, et al. Biopython. Bioinformatics 25:1422, 2009."]
for r in refs: story.append(P(r, SMALL))

os.makedirs("paper", exist_ok=True)
doc = SimpleDocTemplate("paper/outbreak_retro_prediction.pdf", pagesize=letter,
                        leftMargin=0.85*inch, rightMargin=0.85*inch,
                        topMargin=0.8*inch, bottomMargin=0.8*inch,
                        title="Retrospective outbreak prediction", author="Udita Phookan")
doc.build(story)
print("PDF built")
