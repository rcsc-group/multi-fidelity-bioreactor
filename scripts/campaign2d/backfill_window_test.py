"""Backfill window test (lead, 2026-11-04): do the hydrodynamic EDR statistics depend on the scalar release?
For a completed run, compute edr_mean (shear_stress.dat column 9, same nondim scale as hydro_dataset_v2) over
  pre  = the last 3 whole cycles before the tracer release (the edr_mean_t window of the database), and
  post = the first 3 whole cycles after the release.
They must agree within 2 % for the backfill (restart after the release) to be usable as it stands.
Also reports tau_mean (column 5) for the same windows. Read-only.
Usage: uv run python scripts/campaign2d/backfill_window_test.py RUN [RUN ...]
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.postprocess import _load_dat, _t_scales  # noqa: E402

N_CYC = 3
RHO = 1.0e3
TOL = 0.02


def windows(run: str):
    d = ROOT / "runs" / run
    p = json.loads((d / "params.json").read_text())
    arr = _load_dat(d / "shear_stress.dat")
    T_bio, T = _t_scales(p)
    a = p.get("geometry", {}).get("a", 0.25)
    U = a / T_bio
    escale, tscale = RHO * U**3 / a, RHO * U**2
    cyc = arr[:, 1] / T
    c_rel = p.get("t_checkpoint", 0.0) / T + p["n_mix_cycles"]
    c_end = math.floor(c_rel + 1e-6)
    out = {}
    c_last = math.floor(cyc.max() + 1e-6)     # the backfill restarts from the LAST checkpoint: late window
    for name, lo, hi in (("pre", c_end - N_CYC, c_end), ("post", c_end, c_end + N_CYC),
                         ("late", c_last - N_CYC, c_last)):
        if cyc.max() < hi:
            out[name] = None
            continue
        m = (cyc >= lo) & (cyc < hi)
        out[name] = dict(edr=float(arr[m, 9].mean() * escale), tau=float(arr[m, 5].mean() * tscale),
                         window=[lo, hi], n=int(m.sum()))
    return out


def main():
    rel_worst = 0.0
    for run in sys.argv[1:]:
        w = windows(run)
        if w["pre"] is None or w["post"] is None or w["late"] is None:
            print(f"{run}: run too short for the windows")
            continue
        de = w["post"]["edr"] / w["pre"]["edr"] - 1.0
        dt = w["post"]["tau"] / w["pre"]["tau"] - 1.0
        dl = w["late"]["edr"] / w["pre"]["edr"] - 1.0
        rel_worst = max(rel_worst, abs(de), abs(dl))
        print(f"   late window {w['late']['window']}: edr {w['late']['edr']:.5g} ({100 * dl:+.2f} % vs pre)")
        print(f"{run}: edr pre {w['pre']['edr']:.5g} post {w['post']['edr']:.5g} ({100 * de:+.2f} %) | "
              f"tau_mean pre {w['pre']['tau']:.5g} post {w['post']['tau']:.5g} ({100 * dt:+.2f} %) | "
              f"windows {w['pre']['window']} / {w['post']['window']}")
    print(f"worst |edr post/pre - 1| = {100 * rel_worst:.2f} %: {'AGREE' if rel_worst <= TOL else 'DISAGREE'} (limit {100 * TOL:.0f} %)")


if __name__ == "__main__":
    main()
