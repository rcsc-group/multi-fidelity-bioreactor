"""Slide 6 row 3 (round 9): training runs per level, and median core-h per simulated cycle over the runs the chain uses.

Per run: cycles = (tlast - tfirst) / T_run, T_run = _t_scales(params) of that run; core-h/cycle = corehc_all / cycles
(cost_table.tsv 2026-11-04 09:22). Chain sets: 40 L8, 16 L9 (10 cold fig9 theta-7 + 6 warm x8_l9), 13 L10.
For the cold fig9 runs the cold cost is used (their warm estimate is not a measured per-cycle cost).
Output: data/row3_r9.json
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


def per_cycle(name):
    r = ROWS[name]
    if r["tlast"] in ("None", "") or r["tfirst"] in ("None", ""):
        return None
    p = json.load(open(REPO / "runs" / name / "params.json"))
    _, T = _t_scales(p)
    cyc = (float(r["tlast"]) - float(r["tfirst"])) / T
    if cyc <= 0:
        return None
    return float(r["corehc_all"]) / cyc, cyc


sets = {
    "L8": [H1.lf_run(r, t) for r in H1.RPMS for t in H1.LF_THETAS],
    "L9": [H2.l9_run(r, t) for r, t in H2.l9_points()],
    "L10": [H1.hf_run(r, t) for r, t in H1.hf_points()],
}
out = {}
for lev, names in sets.items():
    names = [n for n in names if n in ROWS]
    vals = [per_cycle(n) for n in names]
    keep = [(n, v) for n, v in zip(names, vals) if v is not None]
    pc = [v[0] for _, v in keep]
    cyc = [v[1] for _, v in keep]
    q1, med, q3 = np.percentile(pc, [25, 50, 75])
    out[lev] = dict(n_runs=len(names), n_with_cycles=len(keep), median_core_h_per_cycle=float(med),
                    q1=float(q1), q3=float(q3), median_cycles=float(np.median(cyc)))
    print(f"{lev}: runs {len(names)} (with cycles {len(keep)}), median {med:.2f} core-h/cycle, IQR {q1:.2f}-{q3:.2f}, "
          f"median cycles {np.median(cyc):.1f}")
(TALK / "data" / "row3_r9.json").write_text(json.dumps(out, indent=1))
print("wrote data/row3_r9.json")
