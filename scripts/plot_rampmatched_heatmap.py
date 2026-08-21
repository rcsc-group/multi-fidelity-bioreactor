"""Velocity + shear-stress comparison, ours vs upstream, at a settled
snapshot from the apples-to-apples run (ramp forcing matched exactly,
real cs liquid mask on BOTH sides -- jobs 5083674/5083678, diary.md
2026-08-20). Supersedes the earlier single-instant heatmap, which used
the old mismatched ramp and a contaminated f>0.5-only mask.

Writes two figures (each with colorbars):
  - fields (2x2): |u| and tau, ours vs upstream
  - diff (1x2): |u_ours - u_upstream| / U0, |tau_ours - tau_upstream| / (rho1*U0^2)

Nondimensionalized by U0 (the driver's own "initial rotational
velocity" scale, U0 = w_bio_st*Th_max -- already expressed in the
code's native U_bio-based nondim units, so no extra conversion factor
is needed) rather than by the instantaneous field's own mean: the
mean-based relative error is unstable/not comparable across time when
the field itself is near zero (e.g. during the ramp), whereas a fixed
external scale is stable and physically meaningful throughout.
rho1=1 in the code's nondim units, so the stress scale is just U0^2
(dynamic-pressure convention, consistent with how the code's own
momentum equation is nondimensionalized -- see diary.md 2026-08-20).

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
OUT_PATH_FIELDS = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs/kimetal2024/ours_vs_upstream_study/05_ours_vs_upstream_rampmatched_heatmap.png"
OUT_PATH_RELERR = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs/kimetal2024/ours_vs_upstream_study/06_ours_vs_upstream_rampmatched_relerr.png"

rho_w, mu_w = 1.0e3, 1.0e-3
mu_a = 1.81e-5
L_bio, ANGLE, RPM = 0.25, 7.0, 32.5
Th_max = math.radians(ANGLE)
T_per = 60.0 / RPM
Ly = 0.03575 / L_bio
H_bio = 2. * L_bio * Ly  # [FIX 2026-08-20] was missing the factor of 2 the production driver has had since 2026-08-03 (diary.md) -- Ly is a HALF-height ratio, H_bio must be the FULL bag height
V_bio = L_bio/4*(H_bio + 0.5*L_bio*math.tan(Th_max))
U_bio = V_bio/(H_bio*0.5)/T_per
Re_w = rho_w*U_bio*L_bio/mu_w
mur = mu_a/mu_w
mu1 = 1.0/Re_w
mu2 = mur*mu1
T_bio = L_bio / U_bio
w_bio = 2 * math.pi / T_per
w_bio_st = w_bio * T_bio
U0 = w_bio_st * Th_max  # driver's own characteristic velocity scale; already in code-native nondim (U_bio-based) units
P0 = U0 ** 2  # rho1=1 in code units -> dynamic-pressure stress scale


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
diff_speed[valid_speed] = np.abs(speed_ours[valid_speed] - speed_up[valid_speed]) / U0

diff_tau = np.full((N, N), np.nan)
diff_tau[valid_tau] = np.abs(tau_ours[valid_tau] - tau_up[valid_tau]) / P0

print(f"U0={U0:.4g}  P0={P0:.4g}")
print(f"mean |du|/U0: {np.nanmean(diff_speed):.6g}   mean |dtau|/P0: {np.nanmean(diff_tau):.6g}")
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
TEXT = "#0b0b0b"
CMAP_FIELD = "cividis"
CMAP_TAU = "RdBu_r"
CMAP_ERR = "YlOrRd"  # low end is pale, not black -- magma reads as "no data" at zero (user feedback, diary.md 2026-08-21)

tau_lim = 1.5 * tau_scale


def style_ax(ax):
    ax.set_facecolor(BG)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


# ── Figure 1: |u| and tau, ours vs upstream (2x2, one colorbar per row) ──
fig1, axes1 = plt.subplots(2, 2, figsize=(8, 4.8))
fig1.patch.set_facecolor(BG)

row_specs = [
    ("|u|", [speed_ours, speed_up], CMAP_FIELD, dict(vmin=0, vmax=max(np.nanmax(speed_ours), np.nanmax(speed_up)))),
    ("τ", [tau_ours, tau_up], CMAP_TAU, dict(vmin=-tau_lim, vmax=tau_lim)),
]
for r, (ylabel, (a_ours, a_up), cmap, kw) in enumerate(row_specs):
    im = None
    for c, (field, label) in enumerate([(a_ours, "ours"), (a_up, "upstream")]):
        ax = axes1[r, c]
        style_ax(ax)
        im = ax.imshow(field, origin="lower", cmap=cmap, aspect="equal", **kw)
        if r == 0:
            ax.set_title(label, fontsize=10, color="#52514e", pad=6)
    axes1[r, 0].set_ylabel(ylabel, fontsize=13, color=TEXT, rotation=0, labelpad=20, va="center")
    cbar = fig1.colorbar(im, ax=list(axes1[r, :]), fraction=0.035, pad=0.02)
    cbar.ax.tick_params(labelsize=8, color=TEXT, labelcolor=TEXT)
    cbar.outline.set_visible(False)

fig1.savefig(OUT_PATH_FIELDS, dpi=180, facecolor=fig1.get_facecolor(), bbox_inches="tight")
print(f"Saved to {OUT_PATH_FIELDS}")

# ── Figure 2: nondimensionalized diff of |u| and tau (1x2, each with its own colorbar) ──
fig2, axes2 = plt.subplots(1, 2, figsize=(8, 3.2))
fig2.patch.set_facecolor(BG)

err_specs = [
    ("|Δu| / U0", diff_speed, dict(vmin=0, vmax=np.nanpercentile(diff_speed, 99))),
    ("|Δτ| / (ρU0²)", diff_tau, dict(vmin=0, vmax=np.nanpercentile(diff_tau, 99))),
]
for ax, (label, field, kw) in zip(axes2, err_specs):
    style_ax(ax)
    im = ax.imshow(field, origin="lower", cmap=CMAP_ERR, aspect="equal", **kw)
    ax.set_title(label, fontsize=10, color="#52514e", pad=6)
    cbar = fig2.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.ax.tick_params(labelsize=8, color=TEXT, labelcolor=TEXT)
    cbar.outline.set_visible(False)

fig2.savefig(OUT_PATH_RELERR, dpi=180, facecolor=fig2.get_facecolor(), bbox_inches="tight")
print(f"Saved to {OUT_PATH_RELERR}")
