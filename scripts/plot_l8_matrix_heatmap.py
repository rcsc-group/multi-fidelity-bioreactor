"""4x2 heatmap comparing our own 4 L8-matrix binaries at their final
settled snapshot: |u| and tau (nondimensionalized, code-native units),
one row per binary (diary.md 2026-08-21). One shared colorbar per
column (not per-row) so the 4 rows are directly, visually comparable --
that's the whole point of this figure.

Rows: fresh MPI (3.1), fresh OpenMP (3.3/vanilla), restart-recovered
MPI (3.2, theta=3->7 final state), restart-recovered OpenMP (3.4).
Upstream excluded -- these are specifically "our own binaries" per
the user's request.

Usage:
    uv run python scripts/plot_l8_matrix_heatmap.py
"""
import glob
import math
import re

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

N = 256  # fidelity 8
B_ND = 0.143

RUNS = [
    ("fresh, MPI (3.1)", "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_fresh_mpi"),
    ("fresh, OpenMP (3.3/vanilla)", "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_fresh_openmp"),
    ("restart-recovered, MPI (3.2)", "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_chain_mpi_seg1"),
    ("restart-recovered, OpenMP (3.4)", "/oscar/scratch/eaguerov/tmp/l8_matrix/ours_chain_openmp_seg1"),
]
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/experiments/kimetal2024/ours_vs_upstream_study/08_l8_matrix_mpi_vs_openmp_vs_restart.png"

# Same nondimensionalization as the rest of this session's analysis
# (H_bio corrected 2026-08-20 -- see diary.md, was missing its factor of 2).
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


def latest_time(run_dir):
    files = glob.glob(f"{run_dir}/DataRestart_{N}_*_*.txt")
    pat = re.compile(rf"DataRestart_{N}_([^_]+)_\d+\.txt$")
    times = sorted({float(pat.search(f).group(1)) for f in files if pat.search(f)})
    return times[-1]


def load(run_dir, t):
    files = sorted(glob.glob(f"{run_dir}/DataRestart_{N}_{t:.6g}_*.txt"))
    chunks = [np.loadtxt(f) for f in files]
    return np.vstack(chunks)


def fields(run_dir):
    t = latest_time(run_dir)
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
    return speed, tau, t


# crop to true fluid domain (geometry-derived, same for every row)
yy = (np.arange(N) + 0.5) / N - 0.5
rows_idx = np.where(np.abs(yy) < B_ND)[0]
pad = 3
r0, r1 = max(rows_idx[0] - pad, 0), min(rows_idx[-1] + pad, N)


def crop(a):
    return a[r0:r1, :]


all_speed, all_tau, all_t = [], [], []
for label, run_dir in RUNS:
    s, tau, t = fields(run_dir)
    all_speed.append(crop(s))
    all_tau.append(crop(tau))
    all_t.append(t)
    print(f"{label}: t={t:.4f}  max|u|={np.nanmax(s):.4g}  max|tau|={np.nanmax(np.abs(tau)):.4g}")

speed_vmax = max(np.nanmax(s) for s in all_speed)
tau_vlim = max(np.nanpercentile(np.abs(t)[~np.isnan(t)], 99) for t in all_tau)

BG = "#fcfcfb"
TEXT = "#0b0b0b"
CMAP_FIELD = "cividis"
CMAP_TAU = "RdBu_r"

fig, axes = plt.subplots(4, 2, figsize=(9, 8.5))
fig.patch.set_facecolor(BG)

im_u = im_tau = None
for i, (label, _) in enumerate(RUNS):
    ax_u, ax_tau = axes[i, 0], axes[i, 1]
    for ax in (ax_u, ax_tau):
        ax.set_facecolor(BG)
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    im_u = ax_u.imshow(all_speed[i], origin="lower", cmap=CMAP_FIELD, aspect="equal", vmin=0, vmax=speed_vmax)
    im_tau = ax_tau.imshow(all_tau[i], origin="lower", cmap=CMAP_TAU, aspect="equal", vmin=-tau_vlim, vmax=tau_vlim)
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

fig.savefig(OUT_PATH, dpi=170, facecolor=fig.get_facecolor())
print(f"Saved to {OUT_PATH}")
