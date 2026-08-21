"""4x2 heatmap of nondimensionalized DIFFERENCES between pairs of the L8
matrix's 4 binaries (diary.md 2026-08-21), isolating each factor of the
2x2 design (MPI/OpenMP x fresh/restart) independently, each factor
checked twice for a built-in consistency check:

  row 1: fresh,    MPI vs OpenMP   (same absolute t -- trivial alignment)
  row 2: restart,  MPI vs OpenMP   (same absolute t -- trivial alignment)
  row 3: MPI,      fresh vs restart  (different absolute t -- phase-matched:
                                       the flow locks to absolute t mod T_per,
                                       not "time since start" -- diary.md
                                       2026-08-20 (8) -- so comparing fresh's
                                       last snapshot against the restart
                                       snapshot at the SAME phase, not the
                                       same t, is the fair comparison here)
  row 4: OpenMP,   fresh vs restart  (same phase-matching as row 3)

Nondimensionalized by U0/(rho*U0^2), same convention as
plot_rampmatched_heatmap.py.

Usage:
    uv run python scripts/plot_l8_matrix_relerr_heatmap.py
"""
import glob
import math
import re

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

N = 256
B_ND = 0.143
T_PER_ND = 0.6073  # theta=7, 32.5rpm, corrected H_bio (diary.md 2026-08-20 (6))

FRESH_MPI = "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_fresh_mpi"
FRESH_OPENMP = "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_fresh_openmp"
RESTART_MPI = "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_chain_mpi_seg1"
RESTART_OPENMP = "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_chain_openmp_seg1"
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs/kimetal2024/ours_vs_upstream_study/09_l8_matrix_relerr_heatmap.png"

rho_w, mu_w = 1.0e3, 1.0e-3
mu_a = 1.81e-5
L_bio, ANGLE, RPM = 0.25, 7.0, 32.5
Th_max = math.radians(ANGLE)
T_per = 60.0 / RPM
Ly = 0.03575 / L_bio
H_bio = 2. * L_bio * Ly
V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(Th_max))
U_bio = V_bio / (H_bio * 0.5) / T_per
Re_w = rho_w * U_bio * L_bio / mu_w
mur = mu_a / mu_w
mu1 = 1.0 / Re_w
mu2 = mur * mu1
T_bio = L_bio / U_bio
w_bio = 2 * math.pi / T_per
w_bio_st = w_bio * T_bio
U0 = w_bio_st * Th_max
P0 = U0 ** 2


def mu_of_f(f):
    fc = np.clip(f, 0, 1)
    return 1.0 / (fc * (1.0 / mu1 - 1.0 / mu2) + 1.0 / mu2)


def list_times(run_dir):
    files = glob.glob(f"{run_dir}/DataRestart_{N}_*_*.txt")
    pat = re.compile(rf"DataRestart_{N}_([^_]+)_\d+\.txt$")
    return sorted({float(pat.search(f).group(1)) for f in files if pat.search(f)})


def load(run_dir, t):
    files = sorted(glob.glob(f"{run_dir}/DataRestart_{N}_{t:.6g}_*.txt"))
    return np.vstack([np.loadtxt(f) for f in files])


def fields(run_dir, t):
    data = load(run_dir, t)
    dx = 1.0 / N
    ix = np.clip(np.round((data[:, 0] + 0.5 - dx / 2) / dx).astype(int), 0, N - 1)
    iy = np.clip(np.round((data[:, 1] + 0.5 - dx / 2) / dx).astype(int), 0, N - 1)
    ux = np.full((N, N), np.nan); uy = np.full((N, N), np.nan)
    f = np.full((N, N), np.nan); cs = np.full((N, N), np.nan)
    ux[iy, ix] = data[:, 2]; uy[iy, ix] = data[:, 3]
    f[iy, ix] = data[:, 4]; cs[iy, ix] = data[:, 5]
    du_dy = np.full((N, N), np.nan); dv_dx = np.full((N, N), np.nan)
    du_dy[:, 1:-1] = (ux[:, 2:] - ux[:, :-2]) / (2 * dx)
    dv_dx[1:-1, :] = (uy[2:, :] - uy[:-2, :]) / (2 * dx)
    tau = mu_of_f(f) * (du_dy + dv_dx)
    speed = np.sqrt(ux ** 2 + uy ** 2)
    mask = (f > 0.5) & (cs > 0.5) & ~np.isnan(tau)
    return np.where(mask, speed, np.nan), np.where(mask, tau, np.nan)


def nearest_phase_match(t_ref, candidate_times, settle_after=None):
    # Restrict to the SETTLED tail -- searching the full history risks
    # matching a phase against a still-ramping/transient snapshot (caught
    # empirically: an unrestricted search picked t=11.51, right at
    # t_checkpoint=11.46, mid-ramp -- comparing settled-fresh against
    # barely-restarted is not a fair "restart vs fresh" comparison).
    if settle_after is not None:
        candidate_times = [t for t in candidate_times if t >= settle_after]
    phase_ref = t_ref % T_PER_ND
    best = min(candidate_times, key=lambda t: min(abs((t % T_PER_ND) - phase_ref),
                                                    T_PER_ND - abs((t % T_PER_ND) - phase_ref)))
    return best


# ── row 1/2: same absolute t (trivial) ──
t_fresh = list_times(FRESH_MPI)[-1]
t_fresh_o = min(list_times(FRESH_OPENMP), key=lambda t: abs(t - t_fresh))
t_restart = list_times(RESTART_MPI)[-1]
t_restart_o = min(list_times(RESTART_OPENMP), key=lambda t: abs(t - t_restart))

# ── row 3/4: phase-matched (fresh's last snapshot vs restart's closest-phase
# snapshot, restricted to the settled tail: t_checkpoint + ramp_dur + 3 more
# cycles of margin, so the match can't land mid-ramp) ──
T_CHECKPOINT = 11.46203091341247
RAMP_DUR = 3 * T_PER_ND
SETTLE_AFTER = T_CHECKPOINT + RAMP_DUR + 3 * T_PER_ND
t_restart_phase_for_mpi = nearest_phase_match(t_fresh, list_times(RESTART_MPI), settle_after=SETTLE_AFTER)
t_restart_phase_for_omp = nearest_phase_match(t_fresh_o, list_times(RESTART_OPENMP), settle_after=SETTLE_AFTER)

print(f"row1 (fresh MPI vs OpenMP):   t={t_fresh:.4f} vs t={t_fresh_o:.4f}")
print(f"row2 (restart MPI vs OpenMP): t={t_restart:.4f} vs t={t_restart_o:.4f}")
print(f"row3 (MPI fresh vs restart):  t={t_fresh:.4f} (phase {t_fresh%T_PER_ND:.3f}) vs "
      f"t={t_restart_phase_for_mpi:.4f} (phase {t_restart_phase_for_mpi%T_PER_ND:.3f})")
print(f"row4 (OpenMP fresh vs restart): t={t_fresh_o:.4f} (phase {t_fresh_o%T_PER_ND:.3f}) vs "
      f"t={t_restart_phase_for_omp:.4f} (phase {t_restart_phase_for_omp%T_PER_ND:.3f})")

ROWS = [
    ("fresh: MPI vs OpenMP", FRESH_MPI, t_fresh, FRESH_OPENMP, t_fresh_o),
    ("restart: MPI vs OpenMP", RESTART_MPI, t_restart, RESTART_OPENMP, t_restart_o),
    ("MPI: fresh vs restart\n(phase-matched)", FRESH_MPI, t_fresh, RESTART_MPI, t_restart_phase_for_mpi),
    ("OpenMP: fresh vs restart\n(phase-matched)", FRESH_OPENMP, t_fresh_o, RESTART_OPENMP, t_restart_phase_for_omp),
]

yy = (np.arange(N) + 0.5) / N - 0.5
rows_idx = np.where(np.abs(yy) < B_ND)[0]
pad = 3
r0, r1 = max(rows_idx[0] - pad, 0), min(rows_idx[-1] + pad, N)


def crop(a):
    return a[r0:r1, :]


diff_speed_list, diff_tau_list = [], []
for label, dir_a, t_a, dir_b, t_b in ROWS:
    s_a, tau_a = fields(dir_a, t_a)
    s_b, tau_b = fields(dir_b, t_b)
    valid = ~np.isnan(s_a) & ~np.isnan(s_b)
    diff_speed = np.full((N, N), np.nan)
    diff_speed[valid] = np.abs(s_a[valid] - s_b[valid]) / U0
    valid_tau = ~np.isnan(tau_a) & ~np.isnan(tau_b)
    diff_tau = np.full((N, N), np.nan)
    diff_tau[valid_tau] = np.abs(tau_a[valid_tau] - tau_b[valid_tau]) / P0
    diff_speed_list.append(crop(diff_speed))
    diff_tau_list.append(crop(diff_tau))
    print(f"{label}: mean |du|/U0={np.nanmean(diff_speed):.4g}  mean |dtau|/P0={np.nanmean(diff_tau):.4g}")

speed_vmax = max(np.nanpercentile(d[~np.isnan(d)], 99) for d in diff_speed_list)
tau_vmax = max(np.nanpercentile(d[~np.isnan(d)], 99) for d in diff_tau_list)

BG = "#fcfcfb"
TEXT = "#0b0b0b"
CMAP_ERR = "magma"

fig, axes = plt.subplots(4, 2, figsize=(9, 8.5))
fig.patch.set_facecolor(BG)
im_u = im_tau = None
for i, (label, *_ ) in enumerate(ROWS):
    ax_u, ax_tau = axes[i, 0], axes[i, 1]
    for ax in (ax_u, ax_tau):
        ax.set_facecolor(BG)
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    im_u = ax_u.imshow(diff_speed_list[i], origin="lower", cmap=CMAP_ERR, aspect="equal", vmin=0, vmax=speed_vmax)
    im_tau = ax_tau.imshow(diff_tau_list[i], origin="lower", cmap=CMAP_ERR, aspect="equal", vmin=0, vmax=tau_vmax)
    ax_u.set_ylabel(label, fontsize=8.5, color=TEXT, rotation=0, labelpad=8, ha="right", va="center")

axes[0, 0].set_title("|Δu| / U0", fontsize=11, color=TEXT)
axes[0, 1].set_title("|Δτ| / (ρU0²)", fontsize=11, color=TEXT)
fig.tight_layout(rect=[0.02, 0, 0.94, 1])
cbar_u_ax = fig.add_axes([0.955, 0.55, 0.018, 0.35])
cbar_tau_ax = fig.add_axes([0.955, 0.10, 0.018, 0.35])
cb1 = fig.colorbar(im_u, cax=cbar_u_ax, format=mticker.ScalarFormatter(useMathText=True))
cb1.formatter.set_powerlimits((0, 0)); cb1.update_ticks()
cb2 = fig.colorbar(im_tau, cax=cbar_tau_ax, format=mticker.ScalarFormatter(useMathText=True))
cb2.formatter.set_powerlimits((0, 0)); cb2.update_ticks()
for cb in (cb1, cb2):
    cb.ax.tick_params(labelsize=7, color=TEXT, labelcolor=TEXT)
    cb.outline.set_visible(False)

fig.savefig(OUT_PATH, dpi=170, facecolor=fig.get_facecolor())
print(f"Saved to {OUT_PATH}")
