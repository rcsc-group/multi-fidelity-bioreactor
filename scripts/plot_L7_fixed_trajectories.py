"""Full per-cycle-peak trajectories for the L7 p-fix re-sweep (diary.md
2026-09-07) -- shows the dip-then-rise-to-a-common-plateau shape that a
single settle/never-settle number can't communicate.

Usage:
    uv run python scripts/plot_L7_fixed_trajectories.py
"""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
MANIFEST = json.loads((ROOT / "experiments" / "dtheta_L7_check_v2_fixed_manifest.json").read_text())


def t_per_nd(rpm, theta):
    omega_b = rpm * 2 * math.pi / 60.0
    L_bio, H_bio = 0.25, 2 * 0.03575
    th = math.radians(theta)
    T_per = 2 * math.pi / omega_b
    V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(th))
    U_bio = V_bio / (H_bio * 0.5) / T_per
    T_bio = L_bio / U_bio
    return T_per / T_bio


T_ND = t_per_nd(32.5, 7.0)


def per_cycle_peak(t, tau_mean, t0, n_cycles):
    cycles = (t - t0) / T_ND
    out = np.full(n_cycles, np.nan)
    for c in range(n_cycles):
        mask = (cycles >= c) & (cycles < c + 1)
        if mask.sum() > 0:
            out[c] = tau_mean[mask].max()
    return out


baseline = np.loadtxt(RUNS / "settling_baseline_L7" / "shear_stress.dat", skiprows=1)
n_base = int(baseline[-1, 1] / T_ND)
base_peaks = per_cycle_peak(baseline[:, 1], baseline[:, 5], 0.0, n_base)
target_steady = float(np.nanmean(base_peaks[-5:]))

labels = {"7": r"$\Delta\theta=0$ (trivial)", "6.9": r"$\Delta\theta=0.1°$",
          "4": r"$\Delta\theta=3°$", "2": r"$\Delta\theta=5°$"}
colors = {"7": "#555555", "6.9": "#2166ac", "4": "#b2182b", "2": "#1a9850"}

fig, ax = plt.subplots(figsize=(7, 4.5))
for key in ["7", "6.9", "4", "2"]:
    e = MANIFEST["edges"][key]
    run_id = e["run_id"]
    params = json.loads((RUNS / run_id / "params.json").read_text())
    t_ckpt = params["t_checkpoint"]
    tr = np.loadtxt(RUNS / run_id / "shear_stress.dat", skiprows=1)
    n = int((tr[-1, 1] - t_ckpt) / T_ND)
    peaks = per_cycle_peak(tr[:, 1], tr[:, 5], t_ckpt, n)
    ax.plot(range(n), peaks, color=colors[key], label=labels[key], linewidth=1.3)

ax.axhline(target_steady, color="black", linestyle="--", linewidth=1, label="independent reference")
ax.set_xlabel("Cycles since checkpoint restart")
ax.set_ylabel(r"Peak $\tau_\mathrm{mean}$ per cycle (Pa)")
ax.legend(frameon=False, loc="upper right", fontsize=9)
fig.tight_layout()

out_path = ROOT / "experiments" / "dtheta_monotonicity" / "L7_fixed_trajectories.png"
fig.savefig(out_path, dpi=150)
print(f"Saved {out_path}")
