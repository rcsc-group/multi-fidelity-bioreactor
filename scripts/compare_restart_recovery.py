"""Does the theta=3deg->7deg restart chain converge to the SAME limit
cycle as a genuine fresh theta=7deg cold start? (diary.md 2026-08-20,
user's question). Compares the raw shear_stress.dat/normf.dat time
series over each run's post-ramp tail, phase-folded by t mod T_per_nd,
rather than just the single-number qss/max stats from results.json
(which can hide or exaggerate a phase-offset vs. a genuine amplitude
difference).

Usage:
    uv run python scripts/compare_restart_recovery.py
"""
import math

import numpy as np
import matplotlib.pyplot as plt

RESTART_DIR = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs/c7e9eca7"
BASELINE_DIR = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs/d054ff02"
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/experiments/figures/restart_recovery_comparison.png"

L_bio, ANGLE, RPM = 0.25, 7.0, 32.5
Th_max = math.radians(ANGLE)
T_per = 60.0 / RPM
b_geom = 0.03575
Ly = b_geom / L_bio
H_bio = 2. * L_bio * Ly  # [FIX 2026-08-20] was missing the factor of 2 the production driver has had since 2026-08-03 (diary.md) -- Ly is a HALF-height ratio, H_bio must be the FULL bag height
V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(Th_max))
U_bio = V_bio / (H_bio * 0.5) / T_per
T_bio = L_bio / U_bio
T_per_nd = T_per / T_bio  # at theta=7deg (both runs' actual target condition)


def load_tail(run_dir, n_cycles_tail=8):
    t, tau95, tau98, tau100, taumean = np.loadtxt(f"{run_dir}/shear_stress.dat", skiprows=1,
                                                     usecols=(1, 2, 3, 4, 5), unpack=True)
    _, _, _, _, _, uxrms, _, _, _, uyrms, _, _ = np.loadtxt(f"{run_dir}/normf.dat", skiprows=1,
                                                              usecols=range(2, 14), unpack=True)
    velrms = np.sqrt(uxrms ** 2 + uyrms ** 2)
    t_end = t[-1]
    tail_mask = t > t_end - n_cycles_tail * T_per_nd
    return t[tail_mask], tau95[tail_mask], tau98[tail_mask], velrms[: len(t)][tail_mask]


t_r, tau95_r, tau98_r, vel_r = load_tail(RESTART_DIR)
t_b, tau95_b, tau98_b, vel_b = load_tail(BASELINE_DIR)

phase_r = np.mod(t_r, T_per_nd) / T_per_nd
phase_b = np.mod(t_b, T_per_nd) / T_per_nd

print(f"T_per_nd(theta=7) = {T_per_nd:.4f}")
print(f"restart tail: n={len(t_r)}, t=[{t_r[0]:.2f},{t_r[-1]:.2f}]  tau95 mean={tau95_r.mean():.5g} std={tau95_r.std():.5g}")
print(f"baseline tail: n={len(t_b)}, t=[{t_b[0]:.2f},{t_b[-1]:.2f}]  tau95 mean={tau95_b.mean():.5g} std={tau95_b.std():.5g}")
print(f"restart tail:  velrms mean={vel_r.mean():.5g} std={vel_r.std():.5g}")
print(f"baseline tail: velrms mean={vel_b.mean():.5g} std={vel_b.std():.5g}")

BG = "#fcfcfb"
TEXT = "#0b0b0b"
fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
fig.patch.set_facecolor(BG)
for ax in axes:
    ax.set_facecolor(BG)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

axes[0].scatter(phase_r, tau95_r, s=6, alpha=0.5, color="#2f6f9f", label="restart-recovered (3°→7°)")
axes[0].scatter(phase_b, tau95_b, s=6, alpha=0.5, color="#c1533c", label="fresh (0°→7°)")
axes[0].set_xlabel("phase (t mod T_per / T_per)")
axes[0].set_ylabel("τ₉₅")
axes[0].legend(fontsize=8, frameon=False)

axes[1].scatter(phase_r, vel_r, s=6, alpha=0.5, color="#2f6f9f")
axes[1].scatter(phase_b, vel_b, s=6, alpha=0.5, color="#c1533c")
axes[1].set_xlabel("phase (t mod T_per / T_per)")
axes[1].set_ylabel("vel_rms")

fig.tight_layout()
fig.savefig(OUT_PATH, dpi=160, facecolor=fig.get_facecolor())
print(f"Saved to {OUT_PATH}")
