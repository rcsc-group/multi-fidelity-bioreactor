"""Velocity + shear-stress comparison, ours vs upstream, at a matching
settled L10 snapshot (t=12.7596) -- zero new compute, uses already-run
L10 data (fork_l10_periodic job 5073228, upstream job 4961226).

Usage:
    uv run python scripts/plot_us_vs_upstream_heatmap.py
"""
import glob
import math

import numpy as np
import matplotlib.pyplot as plt

N = 1024
T_SNAPSHOT = 12.7596

OURS_GLOB = f"/oscar/scratch/eaguerov/tmp/fork_l10_periodic/rundir/DataOurs_1024_{T_SNAPSHOT:g}_*.txt"
UPSTREAM_GLOB = f"/oscar/scratch/eaguerov/tmp/upstream_l10/rundir/Data_all/Data_all_1024_{T_SNAPSHOT:g}_*.txt"
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs_site/assets/img/us-vs-upstream-heatmap.png"

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

def fields(glob_pattern):
    data = load_dump(glob_pattern)
    ux, uy, f = to_grid(data)
    dx = 1.0 / N
    du_dy = np.full((N, N), np.nan)
    dv_dx = np.full((N, N), np.nan)
    du_dy[:, 1:-1] = (ux[:, 2:] - ux[:, :-2]) / (2 * dx)
    dv_dx[1:-1, :] = (uy[2:, :] - uy[:-2, :]) / (2 * dx)
    tau = mu_of_f(f) * (du_dy + dv_dx)
    speed = np.sqrt(ux**2 + uy**2)
    mask = f > 0.5
    speed = np.where(mask, speed, np.nan)
    tau = np.where(mask & ~np.isnan(tau), tau, np.nan)
    return speed, tau

speed_ours, tau_ours = fields(OURS_GLOB)
speed_up, tau_up = fields(UPSTREAM_GLOB)

valid_speed = ~np.isnan(speed_ours) & ~np.isnan(speed_up)
valid_tau = ~np.isnan(tau_ours) & ~np.isnan(tau_up)

speed_scale = np.nanmean(speed_up[valid_speed])
tau_scale = np.nanmean(np.abs(tau_up[valid_tau]))

diff_speed = np.full((N, N), np.nan)
diff_speed[valid_speed] = np.abs(speed_ours[valid_speed] - speed_up[valid_speed]) / speed_scale

diff_tau = np.full((N, N), np.nan)
diff_tau[valid_tau] = np.abs(tau_ours[valid_tau] - tau_up[valid_tau]) / tau_scale

# Crop to the actual liquid bounding box (the bag occupies a thin band of
# the full [-0.5,0.5]^2 domain box; most of the array is empty/solid).
rows = np.where(np.any(valid_speed, axis=1))[0]
cols = np.where(np.any(valid_speed, axis=0))[0]
pad = 6
r0, r1 = max(rows[0]-pad, 0), min(rows[-1]+pad, N)
c0, c1 = max(cols[0]-pad, 0), min(cols[-1]+pad, N)

def crop(a):
    return a[r0:r1, c0:c1]

speed_ours, speed_up, diff_speed = crop(speed_ours), crop(speed_up), crop(diff_speed)
tau_ours, tau_up, diff_tau = crop(tau_ours), crop(tau_up), crop(diff_tau)

# ── minimal styling ──
BG = "#fcfcfb"
CMAP_FIELD = "cividis"
CMAP_TAU = "RdBu_r"
CMAP_ERR = "magma"

tau_lim = 1.5 * tau_scale  # tight bulk-scale range -- the 99th pct is
                            # dominated by the known cut-cell singularity and
                            # washes out the bulk pattern entirely; the
                            # outlier just saturates to the colormap extreme.

fig, axes = plt.subplots(2, 3, figsize=(11, 5.2))
fig.patch.set_facecolor(BG)

panels = [
    (speed_ours, "ours", CMAP_FIELD, dict(vmin=0)),
    (speed_up, "upstream", CMAP_FIELD, dict(vmin=0)),
    (diff_speed, "rel. error", CMAP_ERR, dict(vmin=0, vmax=np.nanpercentile(diff_speed, 99))),
    (tau_ours, "ours", CMAP_TAU, dict(vmin=-tau_lim, vmax=tau_lim)),
    (tau_up, "upstream", CMAP_TAU, dict(vmin=-tau_lim, vmax=tau_lim)),
    (diff_tau, "rel. error", CMAP_ERR, dict(vmin=0, vmax=np.nanpercentile(diff_tau, 99))),
]

for ax, (field, label, cmap, kw) in zip(axes.flat, panels):
    ax.set_facecolor(BG)
    ax.imshow(field, origin="lower", cmap=cmap, aspect="equal", **kw)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(label, fontsize=10, color="#52514e", pad=6)

axes[0, 0].set_ylabel("|u|", fontsize=13, color="#0b0b0b", rotation=0, labelpad=20, va="center")
axes[1, 0].set_ylabel("τ", fontsize=13, color="#0b0b0b", rotation=0, labelpad=20, va="center")

fig.tight_layout()
fig.savefig(OUT_PATH, dpi=180, facecolor=fig.get_facecolor())
print(f"Saved to {OUT_PATH}")
print(f"speed_scale={speed_scale:.4g}  tau_scale={tau_scale:.4g}")
print(f"mean rel err speed: {np.nanmean(diff_speed)*100:.1f}%   mean rel err tau: {np.nanmean(diff_tau)*100:.1f}%")
