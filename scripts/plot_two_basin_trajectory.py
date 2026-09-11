"""Static version of the two-basin trajectory evidence (per-cycle peak
tau_mean, cold vs cross-level warm start, normalized to the cold start's
own converged value). Not a video -- a single figure showing the full
60-cycle trajectory at once.

Usage:
    uv run python scripts/plot_two_basin_trajectory.py <warm_run> <cold_run> <out.png> <label>
"""
import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def t_per_nd(rpm=32.5, theta=7.0):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    U = V / (H * 0.5) / T
    return T / (L / U)


def cycle_peak(run):
    d = np.loadtxt(Path("runs") / run / "shear_stress.dat", skiprows=1)
    t, y = d[:, 1], d[:, 5]
    T_ND = t_per_nd()
    c = (t - t[0]) / T_ND
    n = int(c[-1])
    return np.array([y[(c >= k) & (c < k + 1)].max() for k in range(n)])


ap = argparse.ArgumentParser()
ap.add_argument("warm_run")
ap.add_argument("cold_run")
ap.add_argument("out_png")
args = ap.parse_args()

warm = cycle_peak(args.warm_run)
cold = cycle_peak(args.cold_run)
ref = float(np.mean(cold[-15:]))
warm, cold = warm / ref, cold / ref
n = min(len(warm), len(cold))

plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.plot(np.arange(n), cold[:n], color="#3d3d3d", lw=1.8, label="cold start")
ax.plot(np.arange(n), warm[:n], color="#c44536", lw=1.8, label="cross-level warm start")
ax.axhline(1.0, color="#bbbbbb", lw=0.8, zorder=0)
ax.set_xlabel("rocking cycle")
ax.set_ylabel(r"per-cycle peak $\tau_{\rm mean}$ / cold-start converged value")
ax.legend(frameon=False, loc="upper right")
fig.tight_layout()
fig.savefig(args.out_png, dpi=150)
print(f"wrote {args.out_png}")
