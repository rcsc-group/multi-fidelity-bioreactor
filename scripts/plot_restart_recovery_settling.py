"""Raw (non-phase-folded) tau_95(t) for the restart-recovered run vs. the
fresh baseline, plotted since each run's own ramp start -- shows whether
the restart carries a transient overshoot before matching the fresh
run, and how long that takes (diary.md 2026-08-20).

Usage:
    uv run python scripts/plot_restart_recovery_settling.py
"""
import numpy as np
import matplotlib.pyplot as plt

RESTART_DIR = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs/c7e9eca7"
BASELINE_DIR = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs/d054ff02"
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/experiments/figures/restart_recovery_transient_settling.png"


def load(run_dir):
    t, tau95 = np.loadtxt(f"{run_dir}/shear_stress.dat", skiprows=1, usecols=(1, 2), unpack=True)
    return t, tau95


t_r, tau_r = load(RESTART_DIR)
t_b, tau_b = load(BASELINE_DIR)

fig, ax = plt.subplots(figsize=(10, 3.5))
ax.plot(t_r - t_r[0], tau_r, lw=0.6, color="#2f6f9f", label="restart-recovered (3°→7°), t since its own restart")
ax.plot(t_b, tau_b, lw=0.6, color="#c1533c", label="fresh (0°→7°), t since t=0")
ax.set_xlabel("t since ramp start (non-dim)")
ax.set_ylabel("τ₉₅")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT_PATH, dpi=140)
print(f"Saved to {OUT_PATH}")
