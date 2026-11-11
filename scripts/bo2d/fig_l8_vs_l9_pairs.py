"""L9 vs L8 at the SAME site (2026-11-11, to explain "does L8 carry L9 information"). Two panels: kLa_1T_10, tau_95.
Pairs: theta = 7 line (truth snapshots, L8 and L9 at the same rpm, median per site) = filled; 12 off-line pair sites
(pr8_l8 vs bo1_l9, results.json) = hollow. Log-log. No annotations.
  uv run python scripts/bo2d/fig_l8_vs_l9_pairs.py -> figures/bo2d/l8_vs_l9_pairs.png
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402
from scripts.bo2d.pair_l8 import SITES  # noqa: E402

S = ROOT / "scripts/bo2d/state"
KEY = {"kLa_1T_10": ("kLa_1T_10", r"$k_La$ (1/h)"), "tau_95": ("tau_95_qss", r"$\tau_{95}$ (Pa)")}


def line(q):
    d = pd.read_csv(S / f"truth_snapshot_{q}.csv")
    d = d[np.isclose(d.theta, 7) & np.isfinite(d.y) & (d.y > 0)].copy()
    d["r"] = d.rpm.round(1)
    g = d.groupby(["level", "r"]).y.median().unstack(0)
    g = g[[8, 9]].dropna()
    return g[8].to_numpy(), g[9].to_numpy()


def pairs(q):
    k = KEY[q][0] if q == "tau_95" else q
    a, b = [], []
    for x in SITES:
        r8 = json.loads((ROOT / "runs" / f"pr8_l8_rpm{x[0]:.2f}_th{x[1]:.2f}" / "results.json").read_text())
        r9 = json.loads((ROOT / "runs" / f"bo1_l9_rpm{x[0]:.2f}_th{x[1]:.2f}" / "results.json").read_text())
        a.append(r8.get(k, np.nan)); b.append(r9.get(k, np.nan))
    return np.array(a, float), np.array(b, float)


plt.rcParams.update(fs.rcparams())
fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.2))
c8 = fs.level_colour(9)
for ax, q in zip(axes, ("kLa_1T_10", "tau_95")):
    x, y = line(q)
    ax.plot(x, y, "o", ms=7, color=c8, mfc=c8)
    x2, y2 = pairs(q)
    ax.plot(x2, y2, "o", ms=7, color=c8, mfc="w", mew=1.5)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(**fs.GRID_KW)
    ax.set_xlabel(f"L8 {KEY[q][1]}"); ax.set_ylabel(f"L9 {KEY[q][1]}")
from matplotlib.lines import Line2D  # noqa: E402
h = [Line2D([], [], marker="o", ls="none", ms=7, color=c8, mfc=c8), Line2D([], [], marker="o", ls="none", ms=7, color=c8, mfc="w", mew=1.5)]
fig.legend(h, [r"$\theta$ = 7 line", "off the line (12 sites)"], loc="upper center", ncol=2, frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.92))
out = ROOT / "figures/bo2d/l8_vs_l9_pairs.png"
fig.savefig(out, dpi=200)
print("wrote", out)
