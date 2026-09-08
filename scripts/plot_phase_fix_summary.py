"""Quasi-steady offset vs Delta_theta_max, before and after the restart
forcing-phase fix.

Usage:
    uv run python scripts/plot_phase_fix_summary.py
"""
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
OUT = ROOT / "experiments" / "dtheta_monotonicity" / "phase_fix_summary.png"


def t_per_nd(rpm, theta):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    return T / (L / (V / (H * 0.5) / T))


T_ND = t_per_nd(32.5, 7.0)


def late_ratio(run, ref):
    d = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    t, y = d[:, 1], d[:, 5]
    c = (t - t[0]) / T_ND
    p = np.array([y[(c >= k) & (c < k + 1)].max() for k in range(int(c[-1]))])
    return np.mean(p[-20:]) / ref, np.std(p[-20:]) / ref


d = np.loadtxt(RUNS / "kicktest_L7_th7" / "shear_stress.dat", skiprows=1)
c = (d[:, 1] - d[0, 1]) / T_ND
rp = np.array([d[:, 5][(c >= k) & (c < k + 1)].max() for k in range(int(c[-1]))])
REF = float(np.mean(rp[-20:]))

BEFORE = {0.0: "8fd81f04", 0.1: "c2001da9", 3.0: "1ddbba85", 5.0: "3324f34b"}
AFTER = {0.0: "2e69509a", 0.1: "3537850c", 3.0: "8108abc0", 5.0: "fa718f3a"}

fig, ax = plt.subplots(figsize=(6.4, 4.2))
for series, colour, marker, label in [(BEFORE, "#c44536", "o", "before"),
                                      (AFTER, "#2a6f97", "s", "after")]:
    xs = sorted(series)
    ys, es = zip(*(late_ratio(series[x], REF) for x in xs))
    ax.errorbar(xs, ys, yerr=es, color=colour, marker=marker, ms=6, lw=1.3,
                capsize=3, label=label)
ax.axhline(1.0, color="#999999", lw=0.7, zorder=0)
ax.set_xlabel(r"$\Delta\theta_{\max}$ between source and target  [deg]")
ax.set_ylabel(r"quasi-steady $\tau_{\rm mean}$ / cold-start value")
ax.set_title(r"L7, 32.5 rpm, $\theta_{\max}=7^\circ$")
ax.legend(frameon=False, fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(OUT, dpi=160)
print(f"wrote {OUT}")
