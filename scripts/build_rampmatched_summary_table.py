"""Minimal KPI table (percentile of pointwise |u|/tau relative diff),
redone with the apples-to-apples data (ramp-matched, real cs mask on
both sides -- jobs 5083674/5083678). Supersedes the old table, which
used the mismatched ramp and contaminated f>0.5-only mask (diary.md
2026-08-20).

Usage:
    uv run python scripts/build_rampmatched_summary_table.py
"""
import glob
import math

import numpy as np
import matplotlib.pyplot as plt

N = 1024
T_OURS = 12.7447
T_UPSTREAM = 12.7447400981

OURS_GLOB = f"/oscar/scratch/eaguerov/tmp/fork_l10_rampmatch/rundir/DataOurs_1024_{T_OURS:.6g}_*.txt"
UPSTREAM_GLOB = f"/oscar/scratch/eaguerov/tmp/upstream_l10_video/rundir/Data_all/Data_all_1024_{T_UPSTREAM:.12g}_*.txt"
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs/kimetal2024/ours_vs_upstream_study/03_summary_numbers_table.png"

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
speed_up, tau_up = fields(UPSTREAM_GLOB, 7, dict(ux=2, uy=3, f=4, cs=6))

valid_speed = ~np.isnan(speed_ours) & ~np.isnan(speed_up)
valid_tau = ~np.isnan(tau_ours) & ~np.isnan(tau_up)

speed_scale = np.nanmean(speed_up[valid_speed])
tau_scale = np.nanmean(np.abs(tau_up[valid_tau]))

diff_speed = np.abs(speed_ours[valid_speed] - speed_up[valid_speed]) / speed_scale
diff_tau = np.abs(tau_ours[valid_tau] - tau_up[valid_tau]) / tau_scale

rows = []
for p in [50, 90, 99, 100]:
    su = np.percentile(diff_speed, p) * 100
    st = np.percentile(diff_tau, p) * 100
    def fmt(x):
        return f"{x:.2g}%" if x < 1000 else f"{x:.2g}%".replace("e+0", "e")
    rows.append((f"P{p}", fmt(su), fmt(st)))
    print(f"P{p}: |u|={su:.4g}%  tau={st:.4g}%")

BG = "#fcfcfb"
fig, ax = plt.subplots(figsize=(6, 3))
fig.patch.set_facecolor(BG)
ax.set_facecolor(BG)
ax.axis("off")
tbl = ax.table(
    cellText=rows,
    colLabels=["percentile", "|u| pointwise diff", "τ pointwise diff"],
    cellLoc="left",
    colLoc="left",
    loc="center",
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(13)
for (r, c), cell in tbl.get_celld().items():
    cell.set_edgecolor("#d8d6d0" if r > 0 else "none")
    cell.set_facecolor(BG)
    cell.set_text_props(fontweight="bold" if r == 0 else "normal", color="#0b0b0b")
    cell.set_height(0.16)
    if r == 0:
        cell.visible_edges = "B"
        cell.set_edgecolor("#0b0b0b")
    else:
        cell.visible_edges = "B"
fig.tight_layout()
fig.savefig(OUT_PATH, dpi=180, facecolor=fig.get_facecolor())
print(f"Saved to {OUT_PATH}")
