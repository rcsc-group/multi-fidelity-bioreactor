"""Ours-vs-upstream comparison, apples-to-apples run (diary.md 2026-08-19
(5)/(6)): identical ramp forcing on both sides (job 5083674 fork,
5083678 upstream), real `cs` column on BOTH sides (no analytic mask
reconstruction needed anymore), ~225 matched snapshots at 12
frames/rocking-cycle (was 13 sparse snapshots at ~1.5 cycles apart).

Streams frame-by-frame (never holds more than one timestep's fields in
memory) since there are ~225 snapshots x 2 codes x 1024^2 cells. Writes:
  - fields video: |u| and tau, ours vs upstream (2x2, colorbar per row)
  - diff video: |Δu|/U0, |Δτ|/(rho1*U0^2) (1x2, own fixed colorbar each)
    -- nondimensionalized by the driver's own U0 (not by the
    instantaneous field's own mean), so the scale is stable and
    comparable across the whole video, including during the ramp when
    the raw fields are near zero (see plot_rampmatched_heatmap.py for
    the same convention and diary.md 2026-08-20 for the derivation).
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
OUT_VIDEO_RELERR = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs/kimetal2024/ours_vs_upstream_study/07_ours_vs_upstream_rampmatched_relerr_video.mp4"
OUT_CSV = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/experiments/docs/rampmatched_comparison_stats.csv"

rho_w, mu_w = 1.0e3, 1.0e-3
mu_a = 1.81e-5
L_bio, ANGLE, RPM = 0.25, 7.0, 32.5
Th_max = math.radians(ANGLE)
T_per = 60.0 / RPM
b_geom = 0.03575
Ly = b_geom / L_bio
H_bio = 2. * L_bio * Ly  # [FIX 2026-08-20] was missing the factor of 2 the production driver has had since 2026-08-03 (diary.md) -- Ly is a HALF-height ratio, H_bio must be the FULL bag height
V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(Th_max))
U_bio = V_bio / (H_bio * 0.5) / T_per
Re_w = rho_w * U_bio * L_bio / mu_w
mur = mu_a / mu_w
mu1 = 1.0 / Re_w
mu2 = mur * mu1
T_bio = L_bio / U_bio
w_bio = 2 * math.pi / T_per
w_bio_st = w_bio * T_bio
U0 = w_bio_st * Th_max  # driver's own characteristic velocity scale, already code-native nondim
P0 = U0 ** 2  # rho1=1 in code units


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
diff_speed_samples, diff_tau_samples = [], []
for t1 in sample_ts:
    t2 = dict(pairs)[t1]
    s1, tau1, m1 = fields_and_mask(t1, is_upstream=False)
    s2, tau2, m2 = fields_and_mask(t2, is_upstream=True)
    speed_scale_vals.append(np.nanmean(s2[m2]))
    tau_scale_vals.append(np.nanmean(np.abs(tau2[m2])))
    valid = m1 & m2
    diff_speed_samples.append(np.abs(s1[valid] - s2[valid]) / U0)
    diff_tau_samples.append(np.abs(tau1[valid] - tau2[valid]) / P0)
speed_scale = np.mean(speed_scale_vals)
tau_scale = np.mean(tau_scale_vals)
speed_max = 2.0 * speed_scale
tau_lim = 1.5 * tau_scale
diff_speed_max = np.percentile(np.concatenate(diff_speed_samples), 99)
diff_tau_max = np.percentile(np.concatenate(diff_tau_samples), 99)
print(f"speed_scale={speed_scale:.4g} tau_scale={tau_scale:.4g} U0={U0:.4g} P0={P0:.4g}")
print(f"diff_speed_max={diff_speed_max:.4g} diff_tau_max={diff_tau_max:.4g}")

# ── Stream: compute stats + render a frame for every matched pair ──
BG = "#fcfcfb"
TEXT = "#0b0b0b"
CMAP_FIELD = "cividis"
CMAP_TAU = "RdBu_r"
CMAP_ERR = "YlOrRd"  # low end is pale, not black -- magma reads as "no data" at zero (user feedback, diary.md 2026-08-21)
BOX_COLOR = "#0b0b0b"


def style_ax(ax, h, w):
    ax.set_facecolor(BG)
    ax.add_patch(plt.Rectangle((-0.5, -0.5), w - 1, h - 1, fill=False,
                                 edgecolor=BOX_COLOR, linewidth=1.2))
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlim(-2, w + 1)
    ax.set_ylim(-2, h + 1)


tmpdir = Path(tempfile.mkdtemp(prefix="rampmatched_frames_"))
tmpdir_err = Path(tempfile.mkdtemp(prefix="rampmatched_relerr_frames_"))
rows_out = ["t,n_valid,speed_relerr,tau_corr,tau_sign_agree,tau_bulk_mean_relerr,ramp_confounded"]

for i, (t1, t2) in enumerate(pairs):
    s1, tau1, m1 = fields_and_mask(t1, is_upstream=False)
    s2, tau2, m2 = fields_and_mask(t2, is_upstream=True)
    valid = m1 & m2
    n_valid = int(valid.sum())
    diff_speed = np.full((N, N), np.nan)
    diff_tau = np.full((N, N), np.nan)
    if n_valid > 0:
        # nondim by U0/P0 for the plots (stable, comparable across all frames)
        diff_speed[valid] = np.abs(s1[valid] - s2[valid]) / U0
        diff_tau[valid] = np.abs(tau1[valid] - tau2[valid]) / P0
        # relative-error stats (self-referential, for the settling/stability CSV -- unaffected by the plot's normalization choice)
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

    # -- fields frame (2x2, colorbar per row) --
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 3.4))
    fig.patch.set_facecolor(BG)
    row_specs = [
        ("|u|", [s1c, s2c], CMAP_FIELD, dict(vmin=0, vmax=speed_max)),
        ("τ", [tau1c, tau2c], CMAP_TAU, dict(vmin=-tau_lim, vmax=tau_lim)),
    ]
    for r, (ylabel, (f1, f2), cmap, kw) in enumerate(row_specs):
        im = None
        for c, field in enumerate([f1, f2]):
            ax = axes[r, c]
            style_ax(ax, h, w)
            im = ax.imshow(field, origin="lower", cmap=cmap, aspect="equal", **kw)
            if r == 0:
                ax.set_title(["ours", "upstream"][c], fontsize=10, color="#52514e")
        axes[r, 0].set_ylabel(ylabel, fontsize=12, color=TEXT, rotation=0, labelpad=16, va="center")
        cbar = fig.colorbar(im, ax=list(axes[r, :]), fraction=0.035, pad=0.02)
        cbar.ax.tick_params(labelsize=7, color=TEXT, labelcolor=TEXT)
        cbar.outline.set_visible(False)
    frame_path = tmpdir / f"frame_{i:04d}.png"
    fig.savefig(frame_path, dpi=140, facecolor=fig.get_facecolor())
    plt.close(fig)

    # -- relerr frame (1x2, own colorbar each) --
    diff_speed_c, diff_tau_c = crop(diff_speed), crop(diff_tau)
    fige, axese = plt.subplots(1, 2, figsize=(9.6, 2.6))
    fige.patch.set_facecolor(BG)
    err_specs = [
        ("|Δu| / U0", diff_speed_c, dict(vmin=0, vmax=diff_speed_max)),
        ("|Δτ| / (ρU0²)", diff_tau_c, dict(vmin=0, vmax=diff_tau_max)),
    ]
    for ax, (label, field, kw) in zip(axese, err_specs):
        style_ax(ax, h, w)
        im = ax.imshow(field, origin="lower", cmap=CMAP_ERR, aspect="equal", **kw)
        ax.set_title(label, fontsize=10, color="#52514e")
        cbar = fige.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cbar.ax.tick_params(labelsize=7, color=TEXT, labelcolor=TEXT)
        cbar.outline.set_visible(False)
    frame_path_err = tmpdir_err / f"frame_{i:04d}.png"
    fige.savefig(frame_path_err, dpi=140, facecolor=fige.get_facecolor())
    plt.close(fige)

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

subprocess.run([
    "ffmpeg", "-y", "-framerate", "12",
    "-i", str(tmpdir_err / "frame_%04d.png"),
    "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
    "-c:v", "libx264", "-pix_fmt", "yuv420p",
    OUT_VIDEO_RELERR,
], check=True)
print(f"Saved video to {OUT_VIDEO_RELERR}")
