"""README/docs_site hero video: L10 lab-frame velocity-magnitude field,
reconstructed from a real fresh L10 run's per-cell dumps (diary.md
2026-08-28) rather than re-running the native OpenGL video pipeline
(docs_site/how-to/generate-videos.md) at L10, which would need a new
SLURM submission. Rotates the body-frame velocity field by Th(t) --
exact same driver formula as event acceleration() in src/BioReactor.c --
about the domain origin, matching Basilisk's own
quat={0,0,sin(Th/2),cos(Th/2)} lab-frame camera convention (commit
251951c). An earlier version colored by vorticity instead of speed, but
vorticity in this laminar regime is concentrated in thin near-wall
boundary layers and reads as mostly flat/washed-out across a full loop;
velocity magnitude varies richly across the whole bulk flow instead.

RUN_DIR below is a scratch path and not guaranteed to persist (OSCAR
purges scratch after ~30 days untouched, see CLAUDE.md) -- if it's gone,
either point this at a fresh L10 fresh-condition run's dump directory, or
regenerate the hero video properly via the native pipeline instead.
"""
import glob
import math
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import ndimage

N = 1024
RUN_DIR = "/oscar/scratch/eaguerov/tmp/l10_matrix/fresh_mpi"
OUT_PATH = "docs_site/assets/img/hero-rocking-l10-lab.mp4"

L_bio, ANGLE, RPM = 0.25, 7.0, 32.5
b_geom = 0.03575
Th_max = math.radians(ANGLE)
T_per = 60.0 / RPM
Ly = b_geom / L_bio
H_bio = 2.0 * L_bio * Ly
V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(Th_max))
U_bio = V_bio / (H_bio * 0.5) / T_per
T_bio = L_bio / U_bio
w_bio = 2 * math.pi / T_per
w_bio_st = w_bio * T_bio
T_per_st = T_per / T_bio
N_RAMP_CYCLES = 3
ramp_dur = N_RAMP_CYCLES * T_per_st


def theta_body(t):
    x_ss = min(t / ramp_dur, 1.0)
    alpha = 3 * x_ss**2 - 2 * x_ss**3
    return alpha * Th_max * math.sin(w_bio_st * t)


def list_times(run_dir):
    files = glob.glob(f"{run_dir}/DataRestart_{N}_*_*.txt")
    pat = re.compile(rf"DataRestart_{N}_([^_]+)_\d+\.txt$")
    return sorted({float(pat.search(f).group(1)) for f in files if pat.search(f)})


def load(run_dir, t):
    files = sorted(glob.glob(f"{run_dir}/DataRestart_{N}_{t:.6g}_*.txt"))
    return np.vstack([np.loadtxt(f) for f in files])


def build_fields(run_dir, t):
    data = load(run_dir, t)
    dx = 1.0 / N
    ix = np.clip(np.round((data[:, 0] + 0.5 - dx / 2) / dx).astype(int), 0, N - 1)
    iy = np.clip(np.round((data[:, 1] + 0.5 - dx / 2) / dx).astype(int), 0, N - 1)
    f = np.full((N, N), np.nan)
    cs = np.full((N, N), np.nan)
    ux = np.full((N, N), np.nan)
    uy = np.full((N, N), np.nan)
    f[iy, ix] = data[:, 4]
    cs[iy, ix] = data[:, 5]
    ux[iy, ix] = data[:, 2]
    uy[iy, ix] = data[:, 3]
    speed = np.sqrt(ux**2 + uy**2)
    return f, cs, speed


def main():
    times = list_times(RUN_DIR)
    print(f"n_frames = {len(times)}")

    t0_f, t0_cs, _ = build_fields(RUN_DIR, times[0])
    inside0 = t0_cs > 0.5
    rows0 = np.where(inside0.any(axis=1))[0]
    orig_height_px = rows0.max() - rows0.min()
    rotated_height_px = N * math.sin(abs(Th_max)) + orig_height_px * math.cos(abs(Th_max))
    half_px = int(rotated_height_px / 2 * 1.25)
    cy = N // 2

    # Fix a velocity-magnitude color scale from a sample of frames (liquid
    # region only -- speed in air is a separate, much-faster phase).
    sample_idx = np.linspace(0, len(times) - 1, 10).astype(int)
    vmax_samples = []
    for i in sample_idx:
        f, cs, speed = build_fields(RUN_DIR, times[i])
        liq = (cs > 0.5) & (f > 0.5) & ~np.isnan(speed)
        if liq.any():
            vmax_samples.append(np.nanpercentile(speed[liq], 99))
    speed_vmax = max(vmax_samples)
    print(f"speed_vmax (99th pct sample) = {speed_vmax:.3f}")

    cmap = plt.get_cmap("turbo")
    norm = mcolors.Normalize(vmin=0, vmax=speed_vmax)

    tmpdir = Path(tempfile.mkdtemp(prefix="hero_l10_v2_frames_"))
    for fi, t in enumerate(times):
        f, cs, speed = build_fields(RUN_DIR, t)
        inside = cs > 0.5
        liquid = inside & (f > 0.5)
        air = inside & ~liquid

        img = np.ones((N, N, 4))  # RGBA, default opaque white
        img[..., 3] = 0.0  # transparent outside tank by default
        # air: pale, so velocity in the water reads as the main content
        img[air] = [0.92, 0.94, 0.97, 1.0]
        # liquid: colored by velocity magnitude
        liq_rgba = cmap(norm(speed))
        img[liquid] = liq_rgba[liquid]
        img[~inside] = [1.0, 1.0, 1.0, 0.0]

        th_deg = math.degrees(theta_body(t))
        rot = ndimage.rotate(img, th_deg, reshape=False, order=1, mode="constant", cval=0.0)
        rot[..., 3] = np.clip(rot[..., 3], 0, 1)

        lo, hi = max(cy - half_px, 0), min(cy + half_px, N)
        crop = rot[lo:hi, :]

        fig_w, fig_h = 6.0, 6.0 * crop.shape[0] / crop.shape[1]
        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=200)
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")
        ax.imshow(crop, origin="lower", interpolation="bilinear")
        ax.contour(crop[..., 3], levels=[0.5], colors="black", linewidths=1.2)
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.text(0.015, 0.93, f"t = {t*T_bio:.2f} s", transform=ax.transAxes,
                fontsize=13, color="black", va="top")
        fig.tight_layout(pad=0.1)
        fig.savefig(tmpdir / f"frame_{fi:04d}.png", facecolor="white")
        plt.close(fig)
        if fi % 40 == 0:
            print(f"[{fi+1}/{len(times)}]")

    subprocess.run([
        "ffmpeg", "-y", "-framerate", "15",
        "-i", str(tmpdir / "frame_%04d.png"),
        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        "-c:v", "libx264", "-pix_fmt", "yuv420p",
        OUT_PATH,
    ], check=True)
    print(f"Saved to {OUT_PATH}")


if __name__ == "__main__":
    main()
