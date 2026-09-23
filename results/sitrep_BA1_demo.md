# Mutation-risk situation report: Omicron BA.1 (demo)
Generated: 2026-09-23T16:42:22Z by tools/sitrep.py
Risk context: PROMIS pre-outbreak site-risk scores (validated against realized evolution;
see GATES_LOCKED.md / results/gate_results.json).

## Summary
- Mutations assessed: 6
- HIGH pre-outbreak risk-tier sites (top-decile composite): 3
- WATCH tier (top-quartile): 0

## Site-risk annotation (T2)
```
mutation	position	ref_position	composite_z	percentile	tier
S:K417N	417	404	0.24	75%	background
S:E484K	484	470	1.65	95%	HIGH
S:N501Y	501	487	1.60	95%	HIGH
S:D614G	614	600	-0.80	5%	background
S:P681H	681	-	-	-	no-data
S:Q493R	493	479	1.62	95%	HIGH
```

## Target-impact annotation (T3)
```
mutation	gene	position	target_annotations	impact_tier
S:K417N	S	417	RBD;RBD-epitope-class2	HIGH
S:E484K	S	484	RBD;RBD-epitope-class2	HIGH
S:N501Y	S	501	RBD;RBD-epitope-class1	HIGH
S:D614G	S	614	-	low/none-annotated
S:P681H	S	681	furin-cleavage-site	HIGH
S:Q493R	S	493	RBD;RBD-epitope-class1;RBD-epitope-class2	HIGH
```

## Interpretation notes
HIGH/WATCH tiers rank a site against the validated pre-outbreak risk distribution; a HIGH
mutation at a receptor-contact or major-epitope site warrants phenotypic follow-up
(neutralization assays, growth kinetics) and surveillance prioritization. This report is
computational decision support, not a clinical or policy determination.
