#!/bin/bash
# Runs the full fair-ablation grid; output: logs/mf_ablation_fair.log (CSV rows + # comments)
cd ${MFB_REPO:-/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor}
for q in kLa_1T_10 dtmix_0.25 tau_95; do for n in 4 8; do for o in 0 1; do for l in keep drop; do
  uv run python scripts/bo2d/mf_ablation_fair.py $q $n $o $l 2>/dev/null
done; done; done; done
