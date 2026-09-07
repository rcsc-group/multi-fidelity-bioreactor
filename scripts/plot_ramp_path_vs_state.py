"""The 2x2 that separates the ramp code path from the restored state.

Usage:
    uv run python scripts/plot_ramp_path_vs_state.py
"""
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
OUT = ROOT / "experiments" / "dtheta_monotonicity" / "ramp_path_vs_state.png"


def t_per_nd(rpm, theta):
    omega_b = rpm * 2 * math.pi / 60.0
    L_bio, H_bio = 0.25, 2 * 0.03575
    T_per = 2 * math.pi / omega_b
    V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(math.radians(theta)))
    return T_per / (L_bio / (V_bio / (H_bio * 0.5) / T_per))


T_ND = t_per_nd(32.5, 7.0)


def cycle_peaks(run, col=5):
    d = np.loadtxt(RUNS / run / "shear_stress.dat", skiprows=1)
    t, y = d[:, 1], d[:, col]
    cyc = (t - t[0]) / T_ND
    return np.array([y[(cyc >= c) & (cyc < c + 1)].max() for c in range(int(cyc[-1]))])


ref = cycle_peaks("kicktest_L7_th7")
ref_val = float(np.mean(ref[-20:]))

SERIES = [
    ("cold start",                     "kicktest_L7_th7", "#3d3d3d", "-"),
    ("restart, state matched",         "8fd81f04",        "#1f77b4", "-"),
    ("restart, ramp only",             "ba2f3cc1",        "#2ca02c", "--"),
    ("restart, state mismatch only",   "46196ea2",        "#d62728", "-"),
]

fig, ax = plt.subplots(figsize=(7.0, 4.2))
for label, run, c, ls in SERIES:
    p = cycle_peaks(run)
    ax.plot(np.arange(len(p)), p / ref_val, color=c, ls=ls, lw=1.4, label=label)
ax.axhline(1.0, color="#999999", lw=0.7, zorder=0)
ax.set_xlabel("rocking cycle")
ax.set_ylabel(r"per-cycle peak $\tau_{\rm mean}$ / cold-start converged value")
ax.set_title(r"L7, 32.5 rpm, $\theta_{\max}=7^\circ$")
ax.legend(frameon=False, fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(OUT, dpi=160)
print(f"wrote {OUT}")
