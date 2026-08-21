"""4x2 VIDEO comparing our own 4 L8-matrix binaries (diary.md 2026-08-21):
|u| and tau, one row per binary, animated frame-by-frame through each
run's own dump sequence (frame index aligned across rows, not absolute
time -- fresh and restart-recovered runs start at different absolute
t, so "frame i of each row" is the fair comparison, not "same t").
One shared colorbar per column throughout, computed once from a
sample of frames so it stays meaningful across the whole video.

Usage:
    uv run python scripts/render_l8_matrix_video.py
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

N = 256  # fidelity 8
B_ND = 0.143

RUNS = [
    ("fresh, MPI (3.1)", "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_fresh_mpi"),
    ("fresh, OpenMP (3.3)", "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_fresh_openmp"),
    ("restart, MPI (3.2)", "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_chain_mpi_seg1"),
    ("restart, OpenMP (3.4)", "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_chain_openmp_seg1"),
]
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs/kimetal2024/ours_vs_upstream_study/08_l8_matrix_mpi_vs_openmp_vs_restart_video.mp4"

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


def mu_of_f(f):
    fc = np.clip(f, 0, 1)
    return 1.0 / (fc * (1.0 / mu1 - 1.0 / mu2) + 1.0 / mu2)


def list_times(run_dir):
    files = glob.glob(f"{run_dir}/DataRestart_{N}_*_*.txt")
    pat = re.compile(rf"DataRestart_{N}_([^_]+)_\d+\.txt$")
    return sorted({float(pat.search(f).group(1)) for f in files if pat.search(f)})


def load(run_dir, t):
    files = sorted(glob.glob(f"{run_dir}/DataRestart_{N}_{t:.6g}_*.txt"))
    chunks = [np.loadtxt(f) for f in files]
    return np.vstack(chunks)


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
    speed = np.where(mask, speed, np.nan)
    tau = np.where(mask, tau, np.nan)
    return speed, tau


yy = (np.arange(N) + 0.5) / N - 0.5
rows_idx = np.where(np.abs(yy) < B_ND)[0]
pad = 3
r0, r1 = max(rows_idx[0] - pad, 0), min(rows_idx[-1] + pad, N)


def crop(a):
    return a[r0:r1, :]


times_per_run = [list_times(d) for _, d in RUNS]
n_frames = min(len(t) for t in times_per_run)
print(f"n_frames (min across rows) = {n_frames}")

# ── color scale from a handful of sample frames across all rows ──
sample_idx = np.linspace(0, n_frames - 1, 6).astype(int)
speed_vals, tau_vals = [], []
for (_, run_dir), times in zip(RUNS, times_per_run):
    for i in sample_idx:
        s, tau = fields(run_dir, times[i])
        speed_vals.append(np.nanmax(s))
        tau_vals.append(np.nanpercentile(np.abs(tau)[~np.isnan(tau)], 99))
speed_vmax = max(speed_vals)
tau_vlim = max(tau_vals)
print(f"speed_vmax={speed_vmax:.4g} tau_vlim={tau_vlim:.4g}")

BG = "#fcfcfb"
TEXT = "#0b0b0b"
CMAP_FIELD = "cividis"
CMAP_TAU = "RdBu_r"

tmpdir = Path(tempfile.mkdtemp(prefix="l8matrix_frames_"))
for fi in range(n_frames):
    fig, axes = plt.subplots(4, 2, figsize=(9, 8.5))
    fig.patch.set_facecolor(BG)
    im_u = im_tau = None
    for i, ((label, run_dir), times) in enumerate(zip(RUNS, times_per_run)):
        s, tau = fields(run_dir, times[fi])
        s, tau = crop(s), crop(tau)
        ax_u, ax_tau = axes[i, 0], axes[i, 1]
        for ax in (ax_u, ax_tau):
            ax.set_facecolor(BG)
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
        im_u = ax_u.imshow(s, origin="lower", cmap=CMAP_FIELD, aspect="equal", vmin=0, vmax=speed_vmax)
        im_tau = ax_tau.imshow(tau, origin="lower", cmap=CMAP_TAU, aspect="equal", vmin=-tau_vlim, vmax=tau_vlim)
        ax_u.set_ylabel(label, fontsize=9, color=TEXT, rotation=0, labelpad=8, ha="right", va="center")
    axes[0, 0].set_title("|u|", fontsize=12, color=TEXT)
    axes[0, 1].set_title("τ", fontsize=12, color=TEXT)
    fig.tight_layout(rect=[0.02, 0, 0.94, 1])
    cbar_u_ax = fig.add_axes([0.955, 0.55, 0.018, 0.35])
    cbar_tau_ax = fig.add_axes([0.955, 0.10, 0.018, 0.35])
    cb1 = fig.colorbar(im_u, cax=cbar_u_ax)
    cb2 = fig.colorbar(im_tau, cax=cbar_tau_ax, format=mticker.ScalarFormatter(useMathText=True))
    cb2.formatter.set_powerlimits((0, 0))
    cb2.update_ticks()
    for cb in (cb1, cb2):
        cb.ax.tick_params(labelsize=7, color=TEXT, labelcolor=TEXT)
        cb.outline.set_visible(False)
    frame_path = tmpdir / f"frame_{fi:04d}.png"
    fig.savefig(frame_path, dpi=130, facecolor=fig.get_facecolor())
    plt.close(fig)
    if fi % 30 == 0:
        print(f"[{fi+1}/{n_frames}]")

subprocess.run([
    "ffmpeg", "-y", "-framerate", "12",
    "-i", str(tmpdir / "frame_%04d.png"),
    "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
    "-c:v", "libx264", "-pix_fmt", "yuv420p",
    OUT_PATH,
], check=True)
print(f"Saved video to {OUT_PATH}")
