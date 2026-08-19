"""Velocity + shear-stress comparison VIDEO, ours vs upstream, across all
13 matching L10 snapshots (t=0 to 12.76). Zero new compute -- reuses
fork_l10_periodic (job 5073228) and upstream (job 4961226) data.

Usage:
    uv run python scripts/render_us_vs_upstream_video.py
"""
import glob
import math
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

N = 1024
TIMES = [0, 1.0633, 2.1266, 3.1899, 4.2532, 5.3165, 6.3798, 7.4431, 8.5064, 9.5697, 10.633, 11.6963, 12.7596]

OURS_DIR = "/oscar/scratch/eaguerov/tmp/fork_l10_periodic/rundir"
UPSTREAM_DIR = "/oscar/scratch/eaguerov/tmp/upstream_l10/rundir/Data_all"
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs_site/assets/img/us-vs-upstream-comparison.mp4"

rho_w, mu_w = 1.0e3, 1.0e-3
mu_a = 1.81e-5
L_bio, ANGLE, RPM = 0.25, 7.0, 32.5
a_geom, b_geom = 0.25, 0.03575
Ly = b_geom / L_bio
Th_max = math.radians(ANGLE)
T_per = 60.0 / RPM
H_bio = L_bio * Ly
V_bio = L_bio/4*(H_bio + 0.5*L_bio*math.tan(Th_max))
U_bio = V_bio/(H_bio*0.5)/T_per
Re_w = rho_w*U_bio*L_bio/mu_w
mur = mu_a/mu_w
mu1 = 1.0/Re_w
mu2 = mur*mu1

def mu_of_f(f):
    fc = np.clip(f, 0, 1)
    return 1.0/(fc*(1.0/mu1 - 1.0/mu2) + 1.0/mu2)

def load_dump(glob_pattern):
    files = sorted(glob.glob(glob_pattern))
    if not files:
        return None
    chunks = []
    for fp in files:
        try:
            arr = np.loadtxt(fp)
        except ValueError:
            arr = np.loadtxt(fp, skiprows=1)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        chunks.append(arr[:, :5])
    return np.vstack(chunks)

def to_grid(data, n=N):
    dx = 1.0 / n
    ix = np.clip(np.round((data[:, 0] + 0.5 - dx / 2) / dx).astype(int), 0, n - 1)
    iy = np.clip(np.round((data[:, 1] + 0.5 - dx / 2) / dx).astype(int), 0, n - 1)
    grid_ux = np.full((n, n), np.nan)
    grid_uy = np.full((n, n), np.nan)
    grid_f = np.full((n, n), np.nan)
    grid_ux[iy, ix] = data[:, 2]
    grid_uy[iy, ix] = data[:, 3]
    grid_f[iy, ix] = data[:, 4]
    return grid_ux, grid_uy, grid_f

def fields(data):
    ux, uy, f = to_grid(data)
    dx = 1.0 / N
    du_dy = np.full((N, N), np.nan)
    dv_dx = np.full((N, N), np.nan)
    du_dy[:, 1:-1] = (ux[:, 2:] - ux[:, :-2]) / (2 * dx)
    dv_dx[1:-1, :] = (uy[2:, :] - uy[:-2, :]) / (2 * dx)
    tau = mu_of_f(f) * (du_dy + dv_dx)
    mask = f > 0.5
    speed = np.where(mask, np.sqrt(ux**2 + uy**2), np.nan)
    tau = np.where(mask & ~np.isnan(tau), tau, np.nan)
    return speed, tau, mask

# ── Pass 1: load everything, find global crop box + color scales ──
all_data = {}
global_mask = np.zeros((N, N), dtype=bool)
speed_scale_vals, tau_scale_vals = [], []
for t in TIMES:
    ours = load_dump(f"{OURS_DIR}/DataOurs_1024_{t:g}_*.txt")
    up = load_dump(f"{UPSTREAM_DIR}/Data_all_1024_{t:g}_*.txt")
    if ours is None or up is None:
        continue
    s1, tau1, m1 = fields(ours)
    s2, tau2, m2 = fields(up)
    all_data[t] = (s1, tau1, s2, tau2)
    global_mask |= m1 | m2
    if t > 0:
        speed_scale_vals.append(np.nanmean(s2[m2]))
        tau_scale_vals.append(np.nanmean(np.abs(tau2[m2])))

speed_scale = np.mean(speed_scale_vals)
tau_scale = np.mean(tau_scale_vals)
tau_lim = 1.5 * tau_scale
speed_max = 2.0 * speed_scale

rows = np.where(np.any(global_mask, axis=1))[0]
cols = np.where(np.any(global_mask, axis=0))[0]
pad = 6
r0, r1 = max(rows[0]-pad, 0), min(rows[-1]+pad, N)
c0, c1 = max(cols[0]-pad, 0), min(cols[-1]+pad, N)

def crop(a):
    return a[r0:r1, c0:c1]

# ── Pass 2: render one frame per snapshot ──
BG = "#fcfcfb"
CMAP_FIELD = "cividis"
CMAP_TAU = "RdBu_r"
BOX_COLOR = "#0b0b0b"

tmpdir = Path(tempfile.mkdtemp(prefix="us_vs_upstream_frames_"))
frame_paths = []
for i, t in enumerate(TIMES):
    if t not in all_data:
        continue
    s1, tau1, s2, tau2 = all_data[t]
    s1c, s2c, tau1c, tau2c = crop(s1), crop(s2), crop(tau1), crop(tau2)
    h, w = s1c.shape

    fig, axes = plt.subplots(2, 2, figsize=(9, 4.6))
    fig.patch.set_facecolor(BG)
    panels = [
        (s1c, CMAP_FIELD, dict(vmin=0, vmax=speed_max)),
        (s2c, CMAP_FIELD, dict(vmin=0, vmax=speed_max)),
        (tau1c, CMAP_TAU, dict(vmin=-tau_lim, vmax=tau_lim)),
        (tau2c, CMAP_TAU, dict(vmin=-tau_lim, vmax=tau_lim)),
    ]
    for ax, (field, cmap, kw) in zip(axes.flat, panels):
        ax.set_facecolor(BG)
        ax.imshow(field, origin="lower", cmap=cmap, aspect="equal", **kw)
        # Domain bounding box.
        ax.add_patch(plt.Rectangle((-0.5, -0.5), w - 1, h - 1, fill=False,
                                     edgecolor=BOX_COLOR, linewidth=1.5))
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xlim(-2, w+1)
        ax.set_ylim(-2, h+1)

    axes[0, 0].set_title("ours", fontsize=11, color="#52514e")
    axes[0, 1].set_title("upstream", fontsize=11, color="#52514e")
    axes[0, 0].set_ylabel("|u|", fontsize=13, color="#0b0b0b", rotation=0, labelpad=18, va="center")
    axes[1, 0].set_ylabel("τ", fontsize=13, color="#0b0b0b", rotation=0, labelpad=18, va="center")

    fig.tight_layout()
    frame_path = tmpdir / f"frame_{i:03d}.png"
    fig.savefig(frame_path, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)
    frame_paths.append(frame_path)
    print(f"rendered t={t:.4f} -> {frame_path.name}")

# ── Stitch into video (ffmpeg) ──
subprocess.run([
    "ffmpeg", "-y", "-framerate", "2",
    "-i", str(tmpdir / "frame_%03d.png"),
    "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
    "-c:v", "libx264", "-pix_fmt", "yuv420p",
    OUT_PATH,
], check=True)
print(f"\nSaved video to {OUT_PATH}")
