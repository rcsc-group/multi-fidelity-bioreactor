"""Animated version of plot_l8_matrix_relerr_heatmap.py (diary.md
2026-08-21): same 4 pairwise comparisons (MPI/OpenMP x fresh/restart,
each factor isolated twice), animated across fresh's FULL settled tail
(168 frames -- matches the temporal richness of the other L8-matrix
videos, not compressed to one representative cycle).

Frame index is driven by fresh-MPI's own settled-tail time sequence
(its natural cadence, ~14 cycles). Per frame:
  - fresh-OpenMP's partner time: nearest TIME to fresh-MPI's (same
    clock, no phase-matching needed or wanted -- see diary.md 2026-08-21
    (5), independent phase-matching drifted MPI/OpenMP onto adjacent
    cycles and manufactured a fake difference).
  - restart-MPI's partner time: nearest PHASE to fresh-MPI's, within
    restart-MPI's settled tail only (genuinely different clock).
  - restart-OpenMP's partner time: nearest TIME to THAT restart-MPI
    time (same clock as restart-MPI).

Colormap: YlOrRd, not magma -- magma's zero end is black, which reads
as "no data" rather than "measured, confirmed small" for a metric
where zero is the good/expected outcome in half the rows (user's
direct feedback).

Usage:
    uv run python scripts/render_l8_matrix_relerr_video.py
"""
import glob
import math
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

N = 256
B_ND = 0.143
T_PER_ND = 0.6073
N_RAMP_CYCLES = 3
RAMP_DUR = N_RAMP_CYCLES * T_PER_ND
T_CHECKPOINT = 11.46203091341247
SETTLE_FRESH = RAMP_DUR + 3 * T_PER_ND
SETTLE_RESTART = T_CHECKPOINT + RAMP_DUR + 3 * T_PER_ND

FRESH_MPI = "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_fresh_mpi"
FRESH_OPENMP = "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_fresh_openmp"
RESTART_MPI = "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_chain_mpi_seg1"
RESTART_OPENMP = "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_chain_openmp_seg1"
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs/kimetal2024/ours_vs_upstream_study/09_l8_matrix_relerr_video.mp4"

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


def nearest_time(run_dir, settle_after, t_ref):
    times = [t for t in list_times(run_dir) if t >= settle_after]
    return min(times, key=lambda t: abs(t - t_ref))


def nearest_phase(run_dir, settle_after, t_ref):
    times = [t for t in list_times(run_dir) if t >= settle_after]
    phase_ref = t_ref % T_PER_ND
    return min(times, key=lambda t: min(abs((t % T_PER_ND) - phase_ref),
                                          T_PER_ND - abs((t % T_PER_ND) - phase_ref)))


fresh_mpi_times = [t for t in list_times(FRESH_MPI) if t >= SETTLE_FRESH]
n_frames = len(fresh_mpi_times)
print(f"n_frames = {n_frames} (fresh-MPI's full settled tail)")

ROW_LABELS = ["fresh: MPI vs OpenMP", "restart: MPI vs OpenMP",
              "MPI: fresh vs restart", "OpenMP: fresh vs restart"]

yy = (np.arange(N) + 0.5) / N - 0.5
rows_idx = np.where(np.abs(yy) < B_ND)[0]
pad = 3
r0, r1 = max(rows_idx[0] - pad, 0), min(rows_idx[-1] + pad, N)


def crop(a):
    return a[r0:r1, :]


def diff_fields(dir_a, t_a, dir_b, t_b):
    s_a, tau_a = fields(dir_a, t_a)
    s_b, tau_b = fields(dir_b, t_b)
    valid_s = ~np.isnan(s_a) & ~np.isnan(s_b)
    diff_speed = np.full((N, N), np.nan)
    diff_speed[valid_s] = np.abs(s_a[valid_s] - s_b[valid_s]) / U0
    valid_t = ~np.isnan(tau_a) & ~np.isnan(tau_b)
    diff_tau = np.full((N, N), np.nan)
    diff_tau[valid_t] = np.abs(tau_a[valid_t] - tau_b[valid_t]) / P0
    return crop(diff_speed), crop(diff_tau)


def per_frame_diffs(t_fresh_mpi):
    t_fresh_openmp = nearest_time(FRESH_OPENMP, SETTLE_FRESH, t_fresh_mpi)
    t_restart_mpi = nearest_phase(RESTART_MPI, SETTLE_RESTART, t_fresh_mpi)
    t_restart_openmp = nearest_time(RESTART_OPENMP, SETTLE_RESTART, t_restart_mpi)
    return [
        diff_fields(FRESH_MPI, t_fresh_mpi, FRESH_OPENMP, t_fresh_openmp),
        diff_fields(RESTART_MPI, t_restart_mpi, RESTART_OPENMP, t_restart_openmp),
        diff_fields(FRESH_MPI, t_fresh_mpi, RESTART_MPI, t_restart_mpi),
        diff_fields(FRESH_OPENMP, t_fresh_openmp, RESTART_OPENMP, t_restart_openmp),
    ]


# ── color scale from a sample of frames ──
sample_idx = np.linspace(0, n_frames - 1, 8).astype(int)
all_ds, all_dt = [], []
for i in sample_idx:
    for ds, dt in per_frame_diffs(fresh_mpi_times[i]):
        all_ds.append(np.nanpercentile(ds[~np.isnan(ds)], 99))
        all_dt.append(np.nanpercentile(dt[~np.isnan(dt)], 99))
speed_vmax = max(all_ds)
tau_vmax = max(all_dt)
print(f"speed_vmax={speed_vmax:.4g} tau_vmax={tau_vmax:.4g}")

BG = "#fcfcfb"
TEXT = "#0b0b0b"
CMAP_ERR = "YlOrRd"  # low end is pale, not black -- see module docstring

tmpdir = Path(tempfile.mkdtemp(prefix="l8matrix_relerr_frames_"))
for fi, t_fresh_mpi in enumerate(fresh_mpi_times):
    row_diffs = per_frame_diffs(t_fresh_mpi)
    fig, axes = plt.subplots(4, 2, figsize=(9, 8.5))
    fig.patch.set_facecolor(BG)
    im_u = im_tau = None
    for r, label in enumerate(ROW_LABELS):
        ds, dt = row_diffs[r]
        ax_u, ax_tau = axes[r, 0], axes[r, 1]
        for ax in (ax_u, ax_tau):
            ax.set_facecolor(BG)
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
        im_u = ax_u.imshow(ds, origin="lower", cmap=CMAP_ERR, aspect="equal", vmin=0, vmax=speed_vmax)
        im_tau = ax_tau.imshow(dt, origin="lower", cmap=CMAP_ERR, aspect="equal", vmin=0, vmax=tau_vmax)
        ax_u.set_ylabel(label, fontsize=8.5, color=TEXT, rotation=0, labelpad=8, ha="right", va="center")
    axes[0, 0].set_title("|Δu| / U0", fontsize=11, color=TEXT)
    axes[0, 1].set_title("|Δτ| / (ρU0²)", fontsize=11, color=TEXT)
    fig.suptitle(f"t (fresh-MPI clock) = {t_fresh_mpi:.2f}", fontsize=9, color="#52514e", y=0.99)
    fig.tight_layout(rect=[0.02, 0, 0.94, 0.97])
    cbar_u_ax = fig.add_axes([0.955, 0.55, 0.018, 0.35])
    cbar_tau_ax = fig.add_axes([0.955, 0.10, 0.018, 0.35])
    cb1 = fig.colorbar(im_u, cax=cbar_u_ax, format=mticker.ScalarFormatter(useMathText=True))
    cb1.formatter.set_powerlimits((0, 0)); cb1.update_ticks()
    cb2 = fig.colorbar(im_tau, cax=cbar_tau_ax, format=mticker.ScalarFormatter(useMathText=True))
    cb2.formatter.set_powerlimits((0, 0)); cb2.update_ticks()
    for cb in (cb1, cb2):
        cb.ax.tick_params(labelsize=7, color=TEXT, labelcolor=TEXT)
        cb.outline.set_visible(False)
    frame_path = tmpdir / f"frame_{fi:04d}.png"
    fig.savefig(frame_path, dpi=130, facecolor=fig.get_facecolor())
    plt.close(fig)
    if fi % 20 == 0:
        print(f"[{fi+1}/{n_frames}]")

subprocess.run([
    "ffmpeg", "-y", "-framerate", "12",
    "-i", str(tmpdir / "frame_%04d.png"),
    "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
    "-c:v", "libx264", "-pix_fmt", "yuv420p",
    OUT_PATH,
], check=True)
print(f"Saved video to {OUT_PATH}")
