#!/usr/bin/env python3
"""T4 sitrep: mutation-risk situation report (markdown) for a variant.
Usage: python3 tools/sitrep.py --name "ExampleVariant" --case sarscov2 --mutations S:E484K,S:N501Y,...
Combines T2 variant-interpret + T3 target-impact into a situation report."""
import argparse, subprocess, sys, os, datetime

ap = argparse.ArgumentParser()
ap.add_argument("--name", required=True); ap.add_argument("--case", default="sarscov2")
ap.add_argument("--mutations", required=True); ap.add_argument("--out")
a = ap.parse_args()
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
vi = subprocess.run([sys.executable, f"{root}/tools/variant_interpret.py", "--case", a.case,
                     "--mutations", a.mutations], capture_output=True, text=True).stdout
ti = subprocess.run([sys.executable, f"{root}/tools/target_impact.py",
                     "--mutations", a.mutations], capture_output=True, text=True).stdout
n_high = sum(1 for line in vi.splitlines()[1:] if line.endswith("\tHIGH"))
n_watch = sum(1 for line in vi.splitlines()[1:] if line.endswith("\tWATCH"))
muts = a.mutations.split(",")
rep = f"""# Mutation-risk situation report: {a.name}
Generated: {datetime.datetime.utcnow().isoformat(timespec='seconds')}Z by tools/sitrep.py
Risk context: PROMIS pre-outbreak site-risk scores (validated against realized evolution;
see GATES_LOCKED.md / results/gate_results.json).

## Summary
- Mutations assessed: {len(muts)}
- HIGH pre-outbreak risk-tier sites (top-decile composite): {n_high}
- WATCH tier (top-quartile): {n_watch}

## Site-risk annotation (T2)
```
{vi}```

## Target-impact annotation (T3)
```
{ti}```

## Interpretation notes
HIGH/WATCH tiers rank a site against the validated pre-outbreak risk distribution; a HIGH
mutation at a receptor-contact or major-epitope site warrants phenotypic follow-up
(neutralization assays, growth kinetics) and surveillance prioritization. This report is
computational decision support, not a clinical or policy determination.
"""
if a.out: open(a.out, "w").write(rep); print("written", a.out)
else: print(rep)
