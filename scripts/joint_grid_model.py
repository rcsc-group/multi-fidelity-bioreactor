"""Joint multi-level grid-convergence model for mixing time across rpm.

    q_N(x) = q_inf(x) + C(x) h_N^p + e,   h_N = 2^(10-N),   e ~ N(0, s^2)

per chi threshold, fitted to every rpm and level at once. q_inf(x_j), C(x_j)
are free per rpm. p is SHARED, because the convergence order belongs to the scheme, not the
condition. s comes from the residuals, so no noise is assumed. Uncertainty:
the Gauss-Newton covariance s^2 (J^T J)^-1 for every parameter, and 95% = 1.96 sd.
A per-rpm fit with 4 points has no dof to establish p; here p is estimated from
all rpm (10 rpm x 4 levels = 40 points, 21 parameters).

Level sets compared: L6-L9, L7-L9 (drop the coarsest), L6-L8 (cheap only).
Falsifier for "p is common": the residual s is comparable to the level-to-level
changes, i.e. the shared power law does not describe the data.

Usage: uv run python scripts/joint_grid_model.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
RUNS = {9: "fig9_l9_rpm{:g}", 8: "kmix_l8_rpm{:g}", 7: "kmix_l7_rpm{:g}", 6: "kmix_l6_rpm{:g}"}
KIM = pd.read_csv(ROOT / "experiments/kimetal2024/csv_raw/mixing_kla_vs_frequency.csv").set_index("RPM")
SETS = {"L6-L9": (6, 7, 8, 9), "L7-L9": (7, 8, 9), "L6-L8": (6, 7, 8)}


def val(run, key):
    f = ROOT / "runs" / run / "results.json"
    return json.load(open(f)).get(key, math.nan) if f.exists() else math.nan


INVERSE = "--rate" in sys.argv   # fit 1/dtmix (mixing RATE): numerical and physical diffusion add


def data(thr, levels):
    rows = []
    for j, r in enumerate(RPMS):
        for N in levels:
            v = val(RUNS[N].format(r), f"dtmix_{thr:.2f}")
            if np.isfinite(v):
                rows.append((j, 2.0 ** (10 - N), 1.0 / v if INVERSE else v))
    return np.array(rows)


def fit(thr, levels):
    d = data(thr, levels)
    js = sorted(set(d[:, 0].astype(int)))
    # an rpm needs >= 2 levels to separate q_inf from C
    js = [j for j in js if (d[:, 0] == j).sum() >= 2]
    d = d[np.isin(d[:, 0], js)]
    idx = {j: k for k, j in enumerate(js)}
    n = len(js)

    def model(theta):
        qi, C, p = theta[:n], theta[n:2 * n], theta[-1]
        k = np.array([idx[int(j)] for j in d[:, 0]])
        return qi[k] + C[k] * d[:, 1] ** p

    def resid(theta):
        return model(theta) - d[:, 2]

    best = None
    for p0 in (0.5, 1.0, 1.5, 2.0):
        q0 = np.array([d[d[:, 0] == j, 2][np.argmin(d[d[:, 0] == j, 1])] for j in js])
        th0 = np.concatenate([q0, np.full(n, -q0.mean() / 4), [p0]])
        r = least_squares(resid, th0, bounds=(np.r_[np.full(2 * n, -np.inf), 0.05],
                                               np.r_[np.full(2 * n, np.inf), 6.0]))
        if best is None or r.cost < best.cost:
            best = r
    dof = len(d) - (2 * n + 1)
    s2 = 2 * best.cost / dof if dof > 0 else math.nan
    J = best.jac
    cov = s2 * np.linalg.pinv(J.T @ J)
    sd = np.sqrt(np.clip(np.diag(cov), 0, None))
    qi, U = best.x[:n], 1.96 * sd[:n]
    if INVERSE:   # back to a time; delta method for the bound
        qi, U = 1.0 / qi, U / qi ** 2
    out = pd.DataFrame({"rpm": [RPMS[j] for j in js], "q_inf": qi, "U95": U})
    out["kim"] = [float(KIM.loc[r, f"dtmix_strict_{thr:g}"]) for r in out.rpm]
    return out, best.x[-1], 1.96 * sd[-1], math.sqrt(s2), dof, d


def main() -> None:
    for thr in (0.50, 0.75, 0.95):
        print(f"\n===== chi = {thr:.2f}")
        for name, lv in SETS.items():
            out, p, pU, s, dof, d = fit(thr, lv)
            spread = np.median([np.ptp(d[d[:, 0] == j, 2]) for j in set(d[:, 0])])
            ratio = out.q_inf / out.kim
            cover = np.mean(np.abs(out.q_inf - out.kim) <= out.U95)
            print(f"  {name}: p = {p:.2f} +- {pU:.2f}, residual s = {s:.3g} s "
                  f"(median level spread {spread:.3g} s), dof {dof}")
            print(f"     q_inf/Kim {ratio.min():.2f}-{ratio.max():.2f} (median {ratio.median():.2f}), "
                  f"median U95/q_inf {np.median(out.U95 / np.abs(out.q_inf)):.0%}, Kim inside 95%: {cover:.0%}")


if __name__ == "__main__":
    main()
