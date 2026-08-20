"""Phase-binned mean tau_95 profile, restart-recovered vs. fresh
baseline, over their settled tails (diary.md 2026-08-20). Companion
figure to check_restart_recovery_phase_lag.py's printed numbers --
shows the two profiles track each other in phase (no lag needed to
best-match them) with a small, roughly uniform amplitude offset.

Usage:
    uv run python scripts/plot_restart_recovery_phase_profile.py
"""
import numpy as np
import matplotlib.pyplot as plt

RESTART_DIR = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs/c7e9eca7"
BASELINE_DIR = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs/d054ff02"
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/experiments/figures/restart_recovery_phase_profile.png"
T_PER_ND = 0.6073
N_BINS = 16


def load_tail(run_dir, n_cycles_tail=8):
    t, tau95 = np.loadtxt(f"{run_dir}/shear_stress.dat", skiprows=1, usecols=(1, 2), unpack=True)
    t_end = t[-1]
    mask = t > t_end - n_cycles_tail * T_PER_ND
    return t[mask], tau95[mask]


def phase_profile(run_dir):
    t, tau = load_tail(run_dir)
    phase = np.mod(t, T_PER_ND) / T_PER_ND
    bins = np.linspace(0, 1, N_BINS + 1)
    idx = np.clip(np.digitize(phase, bins) - 1, 0, N_BINS - 1)
    mean = np.array([tau[idx == k].mean() for k in range(N_BINS)])
    std = np.array([tau[idx == k].std() for k in range(N_BINS)])
    centers = 0.5 * (bins[:-1] + bins[1:])
    return centers, mean, std


c_r, m_r, s_r = phase_profile(RESTART_DIR)
c_b, m_b, s_b = phase_profile(BASELINE_DIR)

BG = "#fcfcfb"
fig, ax = plt.subplots(figsize=(7, 3.5))
fig.patch.set_facecolor(BG)
ax.set_facecolor(BG)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
ax.errorbar(c_r, m_r, yerr=s_r, fmt="o-", color="#2f6f9f", capsize=2, label="restart-recovered (3°→7°)")
ax.errorbar(c_b, m_b, yerr=s_b, fmt="o-", color="#c1533c", capsize=2, label="fresh (0°→7°)")
ax.set_xlabel("phase (absolute t mod T_per / T_per)")
ax.set_ylabel("τ₉₅ (phase-bin mean ± std)")
ax.legend(fontsize=9, frameon=False)
fig.tight_layout()
fig.savefig(OUT_PATH, dpi=160, facecolor=fig.get_facecolor())
print(f"Saved to {OUT_PATH}")
