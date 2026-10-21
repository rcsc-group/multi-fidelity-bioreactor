"""Settled-state check (pre-registered, diary 2026-10-21): is the L10 training QoI (3-cycle window of the warm-start run) the settled value?

For each L10 point with a 30-cycle continuation (t10_l10_*, or v_l10 at 32.5/7): compare tau_95 (median), tau_mean, EDR of the training window
(scripts/hydro_dataset_v2.qois of the run used in the chain, hf_run) with the same QoI over the LAST three whole cycles of the continuation.
Report per point; median |difference| and max. Restate the claim "within the reference uncertainty" if the median is > 3 % or any point > 8 %.
Usage: uv run python scripts/check_settled_t10.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.hydro_dataset_v2 import COL_EDISS, RHO, qois  # noqa: E402
from scripts.postprocess import _load_dat, _t_scales  # noqa: E402
from scripts.test_h1_mfbml import hf_run  # noqa: E402

POINTS = [(r, 7.0) for r in (17.5, 20, 22.5, 25, 27.5, 30, 32.5)] + [(r, t) for t in (4.0, 9.0) for r in (17.5, 25, 32.5)] + [(r, 5.5) for r in (17.5, 20, 22.5, 27.5, 30, 32.5)]


def late(run, n_cyc=3):
    d = ROOT / "runs" / run
    p = json.load(open(d / "params.json"))
    arr = _load_dat(d / "shear_stress.dat")
    T_bio, T = _t_scales(p)
    U = p.get("geometry", {}).get("a", 0.25) / T_bio
    L = p.get("geometry", {}).get("a", 0.25)
    cyc = arr[:, 1] / T
    c_end = math.floor(cyc.max() + 1e-6)
    m = (cyc >= c_end - n_cyc) & (cyc < c_end)
    return dict(tau_mean_t=float(arr[m, 5].mean() * RHO * U ** 2), tau_95_t=float(np.median(arr[m, 2]) * RHO * U ** 2),
                edr_mean_t=float(arr[m, COL_EDISS].mean() * RHO * U ** 3 / L), c_end=c_end, c_start=float(cyc.min()))


def main():
    rows = []
    print(f"{'rpm/theta':>10} {'cycles':>7} | tau_95 | tau_mean | EDR   (late / training - 1, %)")
    for r, t in POINTS:
        name = "v_l10_rpm32.5_th7" if (r, t) == (32.5, 7.0) else f"t10_l10_rpm{r:g}_th{t:g}"
        run_tr = hf_run(r, t) if t != 5.5 else f"x8_l10_rpm{r:g}_th5.5"
        if not (ROOT / "runs" / name / "shear_stress.dat").exists() or not (ROOT / "runs" / run_tr / "shear_stress.dat").exists():
            continue
        try:
            a, b = late(name), qois(run_tr)
        except Exception as e:  # noqa: BLE001
            print(f"{r:g}/{t:g}: skipped ({e!r:.60})")
            continue
        d = {q: 100 * (a[q] / b[q] - 1) for q in ("tau_95_t", "tau_mean_t", "edr_mean_t")}
        rows.append((r, t, d))
        print(f"{r:5g}/{t:<4g} {a['c_end'] - a['c_start']:7.0f} | {d['tau_95_t']:+6.1f} | {d['tau_mean_t']:+7.1f} | {d['edr_mean_t']:+6.1f}")
    for q in ("tau_95_t", "tau_mean_t", "edr_mean_t"):
        v = np.abs([d[q] for _, _, d in rows])
        print(f"{q}: median |diff| {np.median(v):.1f} %, max {v.max():.1f} % over {len(v)} points")


if __name__ == "__main__":
    main()
