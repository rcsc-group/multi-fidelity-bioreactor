"""Minimal figure: shear-stress percentile sensitivity, upstream's 3 real
settled-state snapshots (already-completed L10 run, zero new compute --
see diary.md 2026-08-18 for where this data comes from and why only
upstream has a multi-point series).

Usage:
    uv run python scripts/plot_tau_percentile_sensitivity.py
"""
import glob

import numpy as np
import matplotlib.pyplot as plt

N = 1024
PERCENTILES = [99.0, 99.9, 100.0]
PERCENTILE_LABELS = ["99th", "99.9th", "100th"]

UPSTREAM_SNAPSHOTS = {
    10.633: "/oscar/scratch/eaguerov/tmp/upstream_l10/rundir/Data_all/Data_all_1024_10.633_*.txt",
    11.6963: "/oscar/scratch/eaguerov/tmp/upstream_l10/rundir/Data_all/Data_all_1024_11.6963_*.txt",
    12.7596: "/oscar/scratch/eaguerov/tmp/upstream_l10/rundir/Data_all/Data_all_1024_12.7596_*.txt",
}
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs_site/assets/img/tau-percentile-sensitivity.png"

# Ordinal blue ramp (dataviz skill palette.md), light->dark for 99th->100th.
COLORS = {"99th": "#86b6ef", "99.9th": "#2a78d6", "100th": "#104281"}
TEXT_SECONDARY = "#52514e"
GRID = "#ececea"


def load_dump(glob_pattern):
    files = sorted(glob.glob(glob_pattern))
    chunks = []
    for fp in files:
        try:
            arr = np.loadtxt(fp)
        except ValueError:
            arr = np.loadtxt(fp, skiprows=1)
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


def tau_percentiles_at(glob_pattern):
    ux, uy, f = to_grid(load_dump(glob_pattern))
    dx = 1.0 / N
    du_dy = np.full((N, N), np.nan)
    dv_dx = np.full((N, N), np.nan)
    du_dy[:, 1:-1] = (ux[:, 2:] - ux[:, :-2]) / (2 * dx)
    dv_dx[1:-1, :] = (uy[2:, :] - uy[:-2, :]) / (2 * dx)
    tau = np.abs(du_dy + dv_dx)
    mask = (f > 0.5) & ~np.isnan(tau)
    return np.percentile(tau[mask], PERCENTILES)


def main():
    t = sorted(UPSTREAM_SNAPSHOTS)
    series = {label: [] for label in PERCENTILE_LABELS}
    for tt in t:
        vals = tau_percentiles_at(UPSTREAM_SNAPSHOTS[tt])
        for label, v in zip(PERCENTILE_LABELS, vals):
            series[label].append(v)

    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    fig.patch.set_facecolor("#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_yscale("log")
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.tick_params(length=0)

    for label in PERCENTILE_LABELS:
        ax.plot(t, series[label], color=COLORS[label], linewidth=2,
                 marker="o", markersize=5, zorder=3)
        ax.annotate(label, (t[-1], series[label][-1]), xytext=(8, 0),
                     textcoords="offset points", va="center",
                     color=COLORS[label], fontsize=10)

    ax.set_xlabel("t", color=TEXT_SECONDARY, fontsize=10)
    ax.set_ylabel("τ", color=TEXT_SECONDARY, fontsize=10, rotation=0, labelpad=10)
    ax.tick_params(colors=TEXT_SECONDARY)
    ax.set_xlim(t[0] - 0.45, t[-1] + 2.25)

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=200, facecolor=fig.get_facecolor())
    print(f"Saved to {OUT_PATH}")


if __name__ == "__main__":
    main()
