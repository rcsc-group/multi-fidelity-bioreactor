"""Velocity + shear-stress comparison, ours vs upstream, at a settled
snapshot from the apples-to-apples run (ramp forcing matched exactly,
real cs liquid mask on BOTH sides -- jobs 5083674/5083678, diary.md
2026-08-20). Supersedes the earlier single-instant heatmap, which used
the old mismatched ramp and a contaminated f>0.5-only mask.

Usage:
    uv run python scripts/plot_rampmatched_heatmap.py
"""
import glob
import math

import numpy as np
import matplotlib.pyplot as plt

N = 1024
B_ND = 0.143
T_OURS = 12.7447   # nearest fork snapshot to the old comparison's t=12.7596
T_UPSTREAM = 12.7447400981  # nearest upstream snapshot to T_OURS (matched within <1e-3)

OURS_GLOB = f"/oscar/scratch/eaguerov/tmp/fork_l10_rampmatch/rundir/DataOurs_1024_{T_OURS:.6g}_*.txt"
UPSTREAM_GLOB_PATTERN = "/oscar/scratch/eaguerov/tmp/upstream_l10_video/rundir/Data_all/Data_all_1024_{t:.12g}_*.txt"
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs/kimetal2024/ours_vs_upstream_study/05_ours_vs_upstream_rampmatched_heatmap.png"

rho_w, mu_w = 1.0e3, 1.0e-3
mu_a = 1.81e-5
L_bio, ANGLE, RPM = 0.25, 7.0, 32.5
Th_max = math.radians(ANGLE)
T_per = 60.0 / RPM
Ly = 0.03575 / L_bio
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


def load(pattern, ncols):
    files = sorted(glob.glob(pattern))
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
    g = {}
    for name, k in cols.items():
        a = np.full((n, n), np.nan)
        a[iy, ix] = data[:, k]
        g[name] = a
    return g


def fields(pattern, ncols, cols):
    data = load(pattern, ncols)
    g = to_grid(data, cols)
    dx = 1.0 / N
    du_dy = np.full((N, N), np.nan)
    dv_dx = np.full((N, N), np.nan)
    du_dy[:, 1:-1] = (g["ux"][:, 2:] - g["ux"][:, :-2]) / (2 * dx)
    dv_dx[1:-1, :] = (g["uy"][2:, :] - g["uy"][:-2, :]) / (2 * dx)
    tau = mu_of_f(g["f"]) * (du_dy + dv_dx)
    speed = np.sqrt(g["ux"] ** 2 + g["uy"] ** 2)
    mask = (g["f"] > 0.5) & (g["cs"] > 0.5)
    speed = np.where(mask, speed, np.nan)
    tau = np.where(mask & ~np.isnan(tau), tau, np.nan)
    return speed, tau


speed_ours, tau_ours = fields(OURS_GLOB, 6, dict(ux=2, uy=3, f=4, cs=5))
speed_up, tau_up = fields(UPSTREAM_GLOB_PATTERN.format(t=T_UPSTREAM), 7, dict(ux=2, uy=3, f=4, cs=6))

valid_speed = ~np.isnan(speed_ours) & ~np.isnan(speed_up)
valid_tau = ~np.isnan(tau_ours) & ~np.isnan(tau_up)

speed_scale = np.nanmean(speed_up[valid_speed])
tau_scale = np.nanmean(np.abs(tau_up[valid_tau]))

diff_speed = np.full((N, N), np.nan)
diff_speed[valid_speed] = np.abs(speed_ours[valid_speed] - speed_up[valid_speed]) / speed_scale

diff_tau = np.full((N, N), np.nan)
diff_tau[valid_tau] = np.abs(tau_ours[valid_tau] - tau_up[valid_tau]) / tau_scale

print(f"mean rel err speed: {np.nanmean(diff_speed)*100:.4f}%   mean rel err tau: {np.nanmean(diff_tau)*100:.4f}%")
print(f"tau corr: {np.corrcoef(tau_ours[valid_tau], tau_up[valid_tau])[0,1]:.4f}")

# Crop to the true fluid domain (geometry-derived, same for both -- see diary.md).
dx = 1.0 / N
yy = (np.arange(N) + 0.5) / N - 0.5
rows = np.where(np.abs(yy) < B_ND)[0]
pad = 6
r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad, N)


def crop(a):
    return a[r0:r1, :]


speed_ours, speed_up, diff_speed = crop(speed_ours), crop(speed_up), crop(diff_speed)
tau_ours, tau_up, diff_tau = crop(tau_ours), crop(tau_up), crop(diff_tau)

BG = "#fcfcfb"
CMAP_FIELD = "cividis"
CMAP_TAU = "RdBu_r"
CMAP_ERR = "magma"

tau_lim = 1.5 * tau_scale

fig, axes = plt.subplots(2, 3, figsize=(11, 4.6))
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
