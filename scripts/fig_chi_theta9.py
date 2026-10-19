"""Mixing curves of the L8 theta-9 grid runs: ln(1 - chi) vs time after release; runs whose 95% mixing time is far from the trend are drawn thick.
Usage: uv run python scripts/fig_chi_theta9.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402
from scripts.outliers_l8 import QS, runs_for, zscores  # noqa: E402
from scripts.test_e8_theta_scatter import GRID, val  # noqa: E402
from scripts.test_s_estimator import curves  # noqa: E402


def main():
    plt.rcParams.update(fs.rcparams())
    runs = runs_for(9.0)
    y = np.array([val(r, "dtmix_0.95") for r in runs])
    z, res = zscores(y)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for r, zi, ri in zip(runs, z, res):
        ts, chi, _ = curves(r)
        m = ts < 220
        far = abs(ri) > 0.35
        ax.plot(ts[m], np.log(np.clip(1 - chi[m], 1e-4, None)), color=fs.level_colour(8), lw=1.6 if far else 0.7, alpha=1.0 if far else 0.45)
    ax.axhline(np.log(0.05), color=fs.GUIDE, lw=1.0)
    ax.set_xlabel("time after release (s)")
    ax.set_ylabel(r"$\ln(1-\chi)$")
    ax.grid(**fs.GRID_KW)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / f"experiments/multifidelity/fig_chi_theta9.{ext}", dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    main()
