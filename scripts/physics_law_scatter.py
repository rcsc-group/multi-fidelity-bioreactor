"""Is the run-to-run scatter of the scalar QoIs explained by the hydrodynamics of the same flow? (zero compute)
L8, 88 grid runs (theta 4, 7, 9; 0.625 rpm grids): ln QoI regressed on a quadratic in (rpm, theta) with and without ln EDR (and ln tau_95, ln tau_mean) from the same run's window.
Leave-one-out RMS ln error (runs within 0.25 rpm at the same theta removed from the training set). QoIs: dtmix_0.50, dtmix_0.95, kLa_1T_10, kLa_1T_25, kLa_1T_50.
Usage: uv run python scripts/physics_law_scatter.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.hydro_dataset_v2 import qois  # noqa: E402

QS = ("dtmix_0.50", "dtmix_0.95", "kLa_1T_10", "kLa_1T_25", "kLa_1T_50")


def load():
    rows = []
    for pat in ("t6_l8_*", "t5_l8_*", "e8_l8_*"):
        for d in sorted((ROOT / "runs").glob(pat)):
            f = d / "results.json"
            if not f.exists():
                continue
            r = json.loads(f.read_text())
            v = [r.get(q) for q in QS]
            try:
                h = qois(d.name)
            except Exception:
                continue
            if h and all(x and x > 0 for x in v):
                p = json.loads((d / "params.json").read_text())
                rows.append([p["omega_b"] * 60 / (2 * math.pi), p["theta_max"][0], *v, h["edr_mean_t"], h["tau_95_t"], h["tau_mean_t"]])
    return np.array(rows)


def design(R, th, hyd=None):
    x, t = (R - 25) / 7.5, (th - 6.5) / 2.5
    cols = [np.ones_like(x), x, t, x * x, x * t, t * t]
    if hyd is not None:
        cols += [hyd] if hyd.ndim == 1 else list(hyd.T)
    return np.column_stack(cols)


def loo(A, y, rpm, th):
    err = []
    for i in range(len(y)):
        tr = ~((np.abs(rpm - rpm[i]) <= 0.25) & (th == th[i]))
        c, *_ = np.linalg.lstsq(A[tr], y[tr], rcond=None)
        err.append(A[i] @ c - y[i])
    return float(np.sqrt(np.mean(np.square(err))))


def main():
    D = load()
    rpm, th = D[:, 0], D[:, 1]
    edr, t95, tm = np.log(D[:, 7]), np.log(D[:, 8]), np.log(D[:, 9])
    edr_c, t95_c, tm_c = (v - v.mean() for v in (edr, t95, tm))
    print(f"L8: {len(D)} runs")
    for j, q in enumerate(QS):
        y = np.log(D[:, 2 + j])
        e0 = loo(design(rpm, th), y, rpm, th)
        e1 = loo(design(rpm, th, edr_c), y, rpm, th)
        e2 = loo(design(rpm, th, np.column_stack([edr_c, t95_c, tm_c])), y, rpm, th)
        print(f"  {q:11s} LOO RMS ln error: (rpm, theta) quadratic {e0:.3f}  + ln EDR {e1:.3f} ({100 * (e1 / e0 - 1):+.0f} %)  + ln EDR, tau_95, tau_mean {e2:.3f} ({100 * (e2 / e0 - 1):+.0f} %)")


if __name__ == "__main__":
    main()
