"""Ours-vs-upstream comparison, apples-to-apples run (diary.md 2026-08-19
(5)/(6)): identical ramp forcing on both sides (job 5083674 fork,
5083678 upstream), real `cs` column on BOTH sides (no analytic mask
reconstruction needed anymore), ~225 matched snapshots at 12
frames/rocking-cycle (was 13 sparse snapshots at ~1.5 cycles apart).

Streams frame-by-frame (never holds more than one timestep's fields in
memory) since there are ~225 snapshots x 2 codes x 1024^2 cells. Writes:
  - a full video (all matched frames)
  - a per-snapshot stats CSV (for checking whether "settled" agreement
    is stable now, and for the phase-lag diagnostic diary.md flagged)

Usage:
    uv run python scripts/analyze_and_render_rampmatched_comparison.py
"""
import glob
import math
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

N = 1024
B_ND = 0.143
T_CHANGE_ST = 11.6132  # upstream's real ramp-completion time (30/T_bio); see diary.md 2026-08-19 (4)

OURS_DIR = "/oscar/scratch/eaguerov/tmp/fork_l10_rampmatch/rundir"
UPSTREAM_DIR = "/oscar/scratch/eaguerov/tmp/upstream_l10_video/rundir/Data_all"
OUT_VIDEO = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs/kimetal2024/ours_vs_upstream_study/04_ours_vs_upstream_rampmatched_video.mp4"
OUT_CSV = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/experiments/docs/rampmatched_comparison_stats.csv"

rho_w, mu_w = 1.0e3, 1.0e-3
mu_a = 1.81e-5
L_bio, ANGLE, RPM = 0.25, 7.0, 32.5
Th_max = math.radians(ANGLE)
T_per = 60.0 / RPM
b_geom = 0.03575
Ly = b_geom / L_bio
H_bio = L_bio * Ly
V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(Th_max))
U_bio = V_bio / (H_bio * 0.5) / T_per
Re_w = rho_w * U_bio * L_bio / mu_w
mur = mu_a / mu_w
mu1 = 1.0 / Re_w
mu2 = mur * mu1


def mu_of_f(f):
    fc = np.clip(f, 0, 1)
    return 1.0 / (fc * (1.0 / mu1 - 1.0 / mu2) + 1.0 / mu2)


def list_times(directory, prefix):
    files = glob.glob(f"{directory}/{prefix}_1024_*_*.txt")
    pat = re.compile(rf"{prefix}_1024_([^_]+)_\d+\.txt$")
    times = sorted({float(pat.search(f).group(1)) for f in files if pat.search(f)})
    return times


def load(glob_pattern, ncols):
    files = sorted(glob.glob(glob_pattern))
    chunks = []
    for fp in files:
        try:
            arr = np.loadtxt(fp)
        except ValueError:
            arr = np.loadtxt(fp, skiprows=1)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        chunks.append(arr[:, :ncols])
    return np.vstack(chunks)


def to_grid(data, cols, n=N):
    dx = 1.0 / n
    ix = np.clip(np.round((data[:, 0] + 0.5 - dx / 2) / dx).astype(int), 0, n - 1)
    iy = np.clip(np.round((data[:, 1] + 0.5 - dx / 2) / dx).astype(int), 0, n - 1)
    grids = {}
    for name, k in cols.items():
        g = np.full((n, n), np.nan)
        g[iy, ix] = data[:, k]
        grids[name] = g
    return grids


def fields_and_mask(t, is_upstream):
    if is_upstream:
        data = load(f"{UPSTREAM_DIR}/Data_all_1024_{t:.12g}_*.txt", 7)
        g = to_grid(data, dict(ux=2, uy=3, f=4, cs=6))
    else:
        data = load(f"{OURS_DIR}/DataOurs_1024_{t:.6g}_*.txt", 6)
        g = to_grid(data, dict(ux=2, uy=3, f=4, cs=5))
    dx = 1.0 / N
    du_dy = np.full((N, N), np.nan)
    dv_dx = np.full((N, N), np.nan)
    du_dy[:, 1:-1] = (g["ux"][:, 2:] - g["ux"][:, :-2]) / (2 * dx)
    dv_dx[1:-1, :] = (g["uy"][2:, :] - g["uy"][:-2, :]) / (2 * dx)
    tau = mu_of_f(g["f"]) * (du_dy + dv_dx)
    speed = np.sqrt(g["ux"] ** 2 + g["uy"] ** 2)
    mask = (g["f"] > 0.5) & (g["cs"] > 0.5) & ~np.isnan(tau)
    return speed, tau, mask


# ── Match snapshot times between the two codes (filenames use different
# float precision, so match nearest within a tight tolerance rather than
# by exact string) ──
ours_times = list_times(OURS_DIR, "DataOurs")
up_times = list_times(UPSTREAM_DIR, "Data_all")
pairs = []
for t1 in ours_times:
    j = np.argmin([abs(t1 - t2) for t2 in up_times])
    if abs(t1 - up_times[j]) < 1e-3:
        pairs.append((t1, up_times[j]))
print(f"{len(ours_times)} ours times, {len(up_times)} upstream times, {len(pairs)} matched pairs")

# ── Crop box from geometry directly (no extra data pass needed): the
# true fluid domain is |y| < B_ND (a_nd doesn't bind -- see diary.md). ──
dx = 1.0 / N
yy = (np.arange(N) + 0.5) / N - 0.5
rows = np.where(np.abs(yy) < B_ND)[0]
pad = 6
r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad, N)
c0, c1 = 0, N  # x never binds; crop only in y for a tight video


def crop(a):
    return a[r0:r1, c0:c1]


# ── Color scale from a handful of settled frames (cheap, avoids a full pass) ──
sample_ts = [t for t, _ in pairs if t >= T_CHANGE_ST][:4]
speed_scale_vals, tau_scale_vals = [], []
for t1 in sample_ts:
    t2 = dict(pairs)[t1]
    s2, tau2, m2 = fields_and_mask(t2, is_upstream=True)
    speed_scale_vals.append(np.nanmean(s2[m2]))
    tau_scale_vals.append(np.nanmean(np.abs(tau2[m2])))
speed_scale = np.mean(speed_scale_vals)
tau_scale = np.mean(tau_scale_vals)
speed_max = 2.0 * speed_scale
tau_lim = 1.5 * tau_scale
print(f"speed_scale={speed_scale:.4g} tau_scale={tau_scale:.4g}")

# ── Stream: compute stats + render a frame for every matched pair ──
BG = "#fcfcfb"
CMAP_FIELD = "cividis"
CMAP_TAU = "RdBu_r"
BOX_COLOR = "#0b0b0b"

tmpdir = Path(tempfile.mkdtemp(prefix="rampmatched_frames_"))
frame_paths = []
rows_out = ["t,n_valid,speed_relerr,tau_corr,tau_sign_agree,tau_bulk_mean_relerr,ramp_confounded"]

for i, (t1, t2) in enumerate(pairs):
    s1, tau1, m1 = fields_and_mask(t1, is_upstream=False)
    s2, tau2, m2 = fields_and_mask(t2, is_upstream=True)
    valid = m1 & m2
    n_valid = int(valid.sum())
    if n_valid > 0:
        speed_relerr = np.mean(np.abs(s1[valid] - s2[valid])) / speed_scale
        a, b = tau1[valid], tau2[valid]
        corr = np.corrcoef(a, b)[0, 1]
        sign_agree = np.mean(np.sign(a) == np.sign(b))
        tau_bulk_mean_relerr = (np.mean(np.abs(a)) - np.mean(np.abs(b))) / tau_scale
    else:
        speed_relerr = corr = sign_agree = tau_bulk_mean_relerr = float("nan")
    confounded = t1 < T_CHANGE_ST
    rows_out.append(f"{t1:.6g},{n_valid},{speed_relerr:.6g},{corr:.6g},{sign_agree:.6g},{tau_bulk_mean_relerr:.6g},{int(confounded)}")

    s1c, s2c, tau1c, tau2c = crop(s1), crop(s2), crop(tau1), crop(tau2)
    h, w = s1c.shape
    fig, axes = plt.subplots(2, 2, figsize=(9, 3.4))
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
        ax.add_patch(plt.Rectangle((-0.5, -0.5), w - 1, h - 1, fill=False,
                                     edgecolor=BOX_COLOR, linewidth=1.2))
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xlim(-2, w + 1)
        ax.set_ylim(-2, h + 1)
    axes[0, 0].set_title("ours", fontsize=10, color="#52514e")
    axes[0, 1].set_title("upstream", fontsize=10, color="#52514e")
    axes[0, 0].set_ylabel("|u|", fontsize=12, color="#0b0b0b", rotation=0, labelpad=16, va="center")
    axes[1, 0].set_ylabel("τ", fontsize=12, color="#0b0b0b", rotation=0, labelpad=16, va="center")
    fig.tight_layout()
    frame_path = tmpdir / f"frame_{i:04d}.png"
    fig.savefig(frame_path, dpi=140, facecolor=fig.get_facecolor())
    plt.close(fig)
    frame_paths.append(frame_path)
    if i % 20 == 0:
        print(f"[{i+1}/{len(pairs)}] t={t1:.4f} speed_relerr={speed_relerr*100:.1f}% tau_corr={corr:+.3f} sign_agree={sign_agree*100:.1f}%")

Path(OUT_CSV).parent.mkdir(parents=True, exist_ok=True)
Path(OUT_CSV).write_text("\n".join(rows_out) + "\n")
print(f"Saved stats CSV to {OUT_CSV}")

subprocess.run([
    "ffmpeg", "-y", "-framerate", "12",
    "-i", str(tmpdir / "frame_%04d.png"),
    "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
    "-c:v", "libx264", "-pix_fmt", "yuv420p",
    OUT_VIDEO,
], check=True)
print(f"Saved video to {OUT_VIDEO}")
