"""Figure: shear-stress percentile sensitivity near the top of the tail,
ours vs. upstream (Kim et al.'s own driver), at a matching settled instant.

Computes tau (naive central-difference stencil, same formula both sides:
mu(f)*(du_dy+dv_dx)) from the raw per-cell velocity dumps captured during
the settled-vs-settled bulk-field comparison (diary.md 2026-08-15/16),
t=12.1466, L10, matching condition -- NOT from shear_stress.dat (which
only logs 95th/98th/100th, not 99th/99.9th, and has no upstream analog).
Shows that a naive stencil blows up at the very top of the percentile
range in BOTH codebases, not just ours -- consistent with the cut-cell
singularity hypothesis being a property of the naive stencil near any
embedded boundary, not a fork-specific bug.

Usage:
    uv run python scripts/plot_tau_percentile_sensitivity.py
"""
import glob

import numpy as np
import matplotlib.pyplot as plt

UPSTREAM_GLOB = "/oscar/scratch/eaguerov/tmp/upstream_l10/rundir/DumpEarly_1024_12.14*_*.txt"
FORK_GLOB = "/oscar/scratch/eaguerov/tmp/fork_l10_coldstart/rundir/DumpEarlyFork_1024_12.14*_*.txt"
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs_site/assets/img/tau-percentile-sensitivity.png"

N = 1024
PERCENTILES = [99.0, 99.9, 100.0]
PERCENTILE_LABELS = ["99th", "99.9th", "100th (max)"]

# Categorical palette (dataviz skill palette.md), fixed order: slot 1 blue (ours),
# slot 2 orange (upstream) -- two distinct datasets, not an ordinal progression.
COLOR_OURS = "#2a78d6"
COLOR_UPSTREAM = "#eb6834"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e3e2dd"


def load_dump(glob_pattern):
    files = sorted(glob.glob(glob_pattern))
    if not files:
        raise SystemExit(f"No files matched: {glob_pattern}")
    chunks = [np.loadtxt(fp)[:, :5] for fp in files]  # x y ux uy f
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


def compute_tau_percentiles(glob_pattern, label):
    print(f"Loading {label}...")
    data = load_dump(glob_pattern)
    ux, uy, f = to_grid(data)
    dx = 1.0 / N

    # Naive central-difference stencil, same formula used everywhere else this
    # session (production event normcal, bio_stress.m): mu=1 (single-phase
    # liquid region only, f>0.5 mask -- mu(f) not reconstructed here since we
    # don't have the two-phase viscosity ratio in this raw dump; consistent
    # relative comparison between the two codebases since both use the same
    # simplification).
    du_dy = np.full((N, N), np.nan)
    dv_dx = np.full((N, N), np.nan)
    du_dy[:, 1:-1] = (ux[:, 2:] - ux[:, :-2]) / (2 * dx)
    dv_dx[1:-1, :] = (uy[2:, :] - uy[:-2, :]) / (2 * dx)
    tau = np.abs(du_dy + dv_dx)

    mask = (f > 0.5) & ~np.isnan(tau)
    tau_liquid = tau[mask]
    print(f"  {label}: {mask.sum()} liquid cells with valid gradient")

    values = np.percentile(tau_liquid, PERCENTILES)
    for p, v in zip(PERCENTILE_LABELS, values):
        print(f"  {label} {p}: {v:.6g}")
    return values


def main():
    ours = compute_tau_percentiles(FORK_GLOB, "ours (fork)")
    upstream = compute_tau_percentiles(UPSTREAM_GLOB, "upstream (Kim et al.)")

    x = np.arange(len(PERCENTILE_LABELS))

    fig, ax = plt.subplots(figsize=(7, 5.5))
    fig.patch.set_facecolor("#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.set_yscale("log")
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)

    ax.plot(x, ours, color=COLOR_OURS, linewidth=2.5, marker="o", markersize=8,
             zorder=3, label="ours (fork)")
    ax.plot(x, upstream, color=COLOR_UPSTREAM, linewidth=2.5, marker="o",
             markersize=8, zorder=3, label="upstream (Kim et al.)")

    ax.annotate("ours", (x[-1], ours[-1]), xytext=(8, 4), textcoords="offset points",
                va="bottom", color=COLOR_OURS, fontsize=11, fontweight="bold")
    ax.annotate("upstream", (x[-1], upstream[-1]), xytext=(8, -10),
                textcoords="offset points", va="top", color=COLOR_UPSTREAM,
                fontsize=11, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(PERCENTILE_LABELS, fontsize=11, color=TEXT_SECONDARY)
    ax.set_xlim(-0.3, len(x) - 1 + 0.55)
    ax.set_ylabel("shear stress |naive stencil|\n(nondimensional, log scale)",
                   color=TEXT_SECONDARY, fontsize=10)
    ax.tick_params(colors=TEXT_SECONDARY)
    ax.set_title(
        "Both codebases blow up at the top of\nthe percentile range",
        color=TEXT_PRIMARY, fontsize=13, fontweight="bold", loc="left", pad=12,
    )
    fig.text(
        0.5, 0.005,
        "L10, θ=7°, 32.5 rpm, settled state (t≈12.15) -- same naive central-difference\n"
        "stencil applied to both codebases' own velocity field",
        ha="center", va="bottom", fontsize=9, color=TEXT_SECONDARY, style="italic",
    )

    fig.tight_layout(rect=(0, 0.055, 1, 1))
    fig.savefig(OUT_PATH, dpi=200, facecolor=fig.get_facecolor())
    print(f"\nSaved to {OUT_PATH}")

    print("\n=== Table version ===")
    print(f"{'percentile':>12} {'ours':>12} {'upstream':>12} {'ours/upstream':>14}")
    for label, o, u in zip(PERCENTILE_LABELS, ours, upstream):
        print(f"{label:>12} {o:>12.6g} {u:>12.6g} {o/u:>14.3f}")


if __name__ == "__main__":
    main()
