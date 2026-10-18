"""Outlier screen of the L8 0.625-rpm grids (theta 4, 7, 9): are the far-off points bad runs or independent randomness?

Per theta and QoI: quadratic trend in log(q) vs rpm on 21 runs; robust z = residual / (1.4826 MAD of residuals); flagged if |z| > 3.
If outliers were bad runs, the same runs would be flagged in several QoIs and their hydrodynamics (tau_95, EDR) would be off too.
Usage: uv run python scripts/outliers_l8.py
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.hydro_dataset_v2 import qois  # noqa: E402
from scripts.test_e8_theta_scatter import GRID, run_id, val  # noqa: E402

QS = ("kLa_1T_10", "kLa_1T_25", "dtmix_0.50", "dtmix_0.95")


def runs_for(th):
    ids = json.loads((ROOT / "experiments/multifidelity/test_t5_roughness.json").read_text())["run_ids"]
    return [ids[f"{r:g}"] for r in GRID] if th == 7.0 else [run_id(r, th) for r in GRID]


def zscores(y):
    x = np.array(GRID)
    ok = np.isfinite(y)
    c = np.polyfit(x[ok], np.log(y[ok]), 2)
    res = np.log(y) - np.polyval(c, x)
    mad = np.nanmedian(np.abs(res[ok] - np.nanmedian(res[ok])))
    return res / (1.4826 * mad), res


def main():
    flags = defaultdict(list)
    out = {}
    for th in (4.0, 7.0, 9.0):
        runs = runs_for(th)
        zs = {}
        for q in QS:
            y = np.array([val(r, q) for r in runs])
            z, res = zscores(y)
            zs[q] = z
            for i in np.where(np.abs(z) > 3)[0]:
                flags[(th, GRID[i])].append((q, float(z[i]), float(np.exp(res[i]))))
        hy = {}
        for h in ("tau_95_t", "edr_mean_t"):
            y = np.array([(qois(r) or {}).get(h, np.nan) for r in runs])
            hy[h], _ = zscores(y)
        out[f"{th:g}"] = dict(z={q: zs[q].tolist() for q in QS}, z_hydro={h: hy[h].tolist() for h in hy})
        print(f"\ntheta {th:g}: flagged runs (|robust z| > 3)")
        for r_i, r in enumerate(GRID):
            fl = flags.get((th, r))
            if fl:
                tau_z, edr_z = hy["tau_95_t"][r_i], hy["edr_mean_t"][r_i]
                print(f"  rpm {r:7.3f} {runs[r_i]:26s} " + "; ".join(f"{q} z {z:+.1f} (x{f:.2f} of trend)" for q, z, f in fl)
                      + f" | hydro z: tau_95 {tau_z:+.1f}, EDR {edr_z:+.1f}")
    n_multi = sum(1 for v in flags.values() if len(v) >= 2)
    print(f"\nruns flagged in >= 2 QoI: {n_multi} of {len(flags)} flagged runs (of 63 runs)")
    (ROOT / "experiments/multifidelity/outliers_l8.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
