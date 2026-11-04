"""Round 9, slide 6 row 3: training runs per level and measured cost per rocking cycle (core-h / cycle).

Per run: cycles simulated = (tlast - tfirst) / T, with T = rocking period in solver time (scripts.postprocess._t_scales).
Cost per cycle = corehc_all / cycles (sacct Elapsed x AllocCPUS summed over all jobs of the run id, cost_table.tsv
2026-11-04 09:22). Level value = median over the runs the chain uses (chain full set: 40 L8, 16 L9, 13 L10).
Output: data/cost_per_cycle_r9.json
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
import scripts.test_h1_mfbml as H1  # noqa: E402
import scripts.test_h2_mfbml as H2  # noqa: E402
from scripts.postprocess import _t_scales  # noqa: E402

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
ROWS = {}
with open("/oscar/scratch/eaguerov/tmp/cost_audit/cost_table.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        ROWS[r["run"]] = r
T_PARAMS = json.load(open(REPO / "runs" / "kmix_l8_rpm32.5_th7" / "params.json")) if (REPO / "runs" / "kmix_l8_rpm32.5_th7").exists() else None
_, T = _t_scales(T_PARAMS) if T_PARAMS else (None, 0.607329)
print("T (solver time per cycle) =", T)

sets = {
    "L8": [H1.lf_run(r, t) for r, t in [(r, t) for r in H1.RPMS for t in H1.LF_THETAS]],
    "L9": [H2.l9_run(r, t) for r, t in H2.l9_points()],
    "L10": [H1.hf_run(r, t) for r, t in H1.hf_points()],
}
out = {}
for lev, names in sets.items():
    names = [n for n in names if n in ROWS]
    vals = []
    for n in names:
        r = ROWS[n]
        cyc = (float(r["tlast"]) - float(r["tfirst"])) / T
        cost = float(r["corehc_all"])
        if cyc > 0:
            vals.append((n, cost, cyc, cost / cyc))
    med = float(np.median([v[3] for v in vals]))
    out[lev] = dict(n_runs=len(names), n_with_cycles=len(vals), median_core_h_per_cycle=med,
                    min=float(min(v[3] for v in vals)), max=float(max(v[3] for v in vals)),
                    median_cycles=float(np.median([v[2] for v in vals])))
    print(f"{lev}: runs {len(names)}, with cycles {len(vals)}, median core-h/cycle {med:.2f} "
          f"(range {out[lev]['min']:.2f} to {out[lev]['max']:.2f}), median cycles {out[lev]['median_cycles']:.1f}")
(TALK / "data" / "cost_per_cycle_r9.json").write_text(json.dumps(out, indent=1))
print("wrote data/cost_per_cycle_r9.json")
