"""Test W (pre-registered, diary 2026-10-18): width of the 22.5 rpm anomaly at L8. New runs w_l8_rpm{21.75, 22, 22.25, 22.75, 23, 23.25}_th{7, 4}.
Robust z of each new run against the quadratic trend and MAD of the existing 21 grid runs. LOCALIZED REGIME if theta 7 EDR z >= 5 in >= 2 of the 6 new runs;
SINGLE-POINT if all 6 new theta-7 runs have |EDR z| <= 3; else MIXED.
Usage: uv run python scripts/test_w_width.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.hydro_dataset_v2 import qois  # noqa: E402
from scripts.outliers_l8 import runs_for  # noqa: E402
from scripts.test_e8_theta_scatter import GRID, val  # noqa: E402

NEW = (21.75, 22.0, 22.25, 22.75, 23.0, 23.25)


def getter(q):
    return (lambda r: (qois(r) or {}).get(q, np.nan)) if q.endswith("_t") else (lambda r: val(r, q))


def main():
    res = {}
    for th in (7.0, 4.0):
        runs = runs_for(th)
        print(f"\ntheta {th:g}: robust z of the new runs (existing grid: 21 runs)")
        print("  rpm      " + "  ".join(f"{q:>10s}" for q in ("edr_mean_t", "tau_95_t", "dtmix_0.50", "kLa_1T_10", "dtmix_0.95")))
        zs = {}
        for q in ("edr_mean_t", "tau_95_t", "dtmix_0.50", "kLa_1T_10", "dtmix_0.95"):
            g = getter(q)
            y = np.array([g(r) for r in runs])
            c = np.polyfit(GRID, np.log(y), 2)
            r0 = np.log(y) - np.polyval(c, GRID)
            mad = 1.4826 * np.median(np.abs(r0 - np.median(r0)))
            zs[q] = [(np.log(g(f"w_l8_rpm{r:g}_th{th:g}")) - np.polyval(c, r)) / mad for r in NEW]
        for i, r in enumerate(NEW):
            print(f"  {r:6.2f}   " + "  ".join(f"{zs[q][i]:+10.1f}" for q in zs))
        res[th] = zs
    z7 = np.array(res[7.0]["edr_mean_t"])
    verdict = "LOCALIZED REGIME" if (z7 >= 5).sum() >= 2 else "SINGLE-POINT" if np.all(np.abs(z7) <= 3) else "MIXED"
    print("\nW verdict (theta 7 EDR):", verdict, " z =", np.round(z7, 1))


if __name__ == "__main__":
    main()
