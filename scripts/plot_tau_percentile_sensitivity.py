"""Figure: shear-stress percentile sensitivity over time, ours vs. upstream,
at L10 -- using ONLY already-computed data from the two L10 jobs we already
ran (upstream job 4877551/4961226, ours job 4961227), zero new compute.

Upstream's own `OUT_FILES` mechanism (unrelated to our comparison -- just a
pre-existing feature of Kim et al.'s driver, `dt_file~=1.06`) dumped periodic
full-field snapshots throughout its ENTIRE L10 run "for free". Three of
those land after upstream's own ~16.25-cycle ramp (t=10.63, 11.70, 12.76),
giving a real (if coarse) 3-point settled time series at zero extra cost.

Our own fork's L10 run was never set up with an equivalent periodic dump --
it only has the single settled snapshot (t=12.15) used for the original
bulk-field comparison (diary.md 2026-08-15/16). Extending it would cost the
same per-step price as any other L10 continuation (~hours), so this figure
shows that single point honestly rather than pretending a matching series
exists. See diary.md 2026-08-18 for the full reasoning.

Usage:
    uv run python scripts/plot_tau_percentile_sensitivity.py
"""
import glob

import numpy as np
import matplotlib.pyplot as plt

N = 1024
PERCENTILES = [99.0, 99.9, 100.0]
PERCENTILE_LABELS = ["99th", "99.9th", "100th (max)"]
PERCENTILE_MARKER = {"99th": "^", "99.9th": "s", "100th (max)": "o"}

COLOR_OURS = "#2a78d6"
COLOR_UPSTREAM = "#eb6834"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e3e2dd"

# Upstream: 3 settled full-field dumps already on disk (past its own ramp,
# t_change_st~=9.87), from OUT_FILES -- x y ux uy f columns 0-4 of 12.
UPSTREAM_SNAPSHOTS = {
    10.633: "/oscar/scratch/eaguerov/tmp/upstream_l10/rundir/Data_all/Data_all_1024_10.633_*.txt",
    11.6963: "/oscar/scratch/eaguerov/tmp/upstream_l10/rundir/Data_all/Data_all_1024_11.6963_*.txt",
    12.7596: "/oscar/scratch/eaguerov/tmp/upstream_l10/rundir/Data_all/Data_all_1024_12.7596_*.txt",
}
# Ours: the single settled snapshot from the bulk-field comparison.
OURS_SNAPSHOT = {
    12.1465832326: "/oscar/scratch/eaguerov/tmp/fork_l10_coldstart/rundir/DumpEarlyFork_1024_12.1465832326_*.txt",
}
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs_site/assets/img/tau-percentile-sensitivity.png"


def load_dump(glob_pattern):
    files = sorted(glob.glob(glob_pattern))
    if not files:
        raise SystemExit(f"No files matched: {glob_pattern}")
    # Data_all_*.txt (upstream OUT_FILES) has a header row; DumpEarlyFork_*.txt
    # (our single snapshot) does not -- skip a possible header transparently.
    chunks = []
    for fp in files:
        try:
            arr = np.loadtxt(fp)
        except ValueError:
            arr = np.loadtxt(fp, skiprows=1)
        chunks.append(arr[:, :5])  # x y ux uy f
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
    data = load_dump(glob_pattern)
    ux, uy, f = to_grid(data)
    dx = 1.0 / N
    du_dy = np.full((N, N), np.nan)
    dv_dx = np.full((N, N), np.nan)
    du_dy[:, 1:-1] = (ux[:, 2:] - ux[:, :-2]) / (2 * dx)
    dv_dx[1:-1, :] = (uy[2:, :] - uy[:-2, :]) / (2 * dx)
    tau = np.abs(du_dy + dv_dx)
    mask = (f > 0.5) & ~np.isnan(tau)
    return np.percentile(tau[mask], PERCENTILES)


def main():
    print("Computing upstream percentiles at each settled snapshot...")
    t_up = sorted(UPSTREAM_SNAPSHOTS)
    pctl_up = {label: [] for label in PERCENTILE_LABELS}
    for t in t_up:
        vals = tau_percentiles_at(UPSTREAM_SNAPSHOTS[t])
        for label, v in zip(PERCENTILE_LABELS, vals):
            pctl_up[label].append(v)
        print(f"  t={t}: {dict(zip(PERCENTILE_LABELS, vals))}")

    print("Computing ours percentile at its single settled snapshot...")
    t_ours = sorted(OURS_SNAPSHOT)
    pctl_ours = {label: [] for label in PERCENTILE_LABELS}
    for t in t_ours:
        vals = tau_percentiles_at(OURS_SNAPSHOT[t])
        for label, v in zip(PERCENTILE_LABELS, vals):
            pctl_ours[label].append(v)
        print(f"  t={t}: {dict(zip(PERCENTILE_LABELS, vals))}")

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor("#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.set_yscale("log")
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)

    for label in PERCENTILE_LABELS:
        ax.plot(t_up, pctl_up[label], color=COLOR_UPSTREAM, linewidth=2,
                 marker=PERCENTILE_MARKER[label], markersize=8, zorder=3)
        ax.annotate(label, (t_up[-1], pctl_up[label][-1]), xytext=(8, 0),
                     textcoords="offset points", va="center", color=COLOR_UPSTREAM,
                     fontsize=9)
        # Ours: a single point per percentile -- no line, honestly reflecting
        # that only one snapshot exists (see module docstring).
        ax.plot(t_ours, pctl_ours[label], color=COLOR_OURS, linestyle="none",
                 marker=PERCENTILE_MARKER[label], markersize=11,
                 markeredgecolor="white", markeredgewidth=1, zorder=4)

    ax.annotate("ours (single L10 snapshot)", (t_ours[0], pctl_ours["100th (max)"][0]),
                 xytext=(10, 10), textcoords="offset points", va="bottom",
                 color=COLOR_OURS, fontsize=10, fontweight="bold")
    ax.annotate("upstream (3 free snapshots\nfrom its own periodic dump)",
                 (t_up[0], pctl_up["100th (max)"][0]), xytext=(-10, 12),
                 textcoords="offset points", va="bottom", ha="right",
                 color=COLOR_UPSTREAM, fontsize=10, fontweight="bold")

    from matplotlib.lines import Line2D
    marker_legend = [
        Line2D([0], [0], color=TEXT_SECONDARY, marker=PERCENTILE_MARKER[l],
               linestyle="none", markersize=8, label=l)
        for l in PERCENTILE_LABELS
    ]
    leg = ax.legend(handles=marker_legend, loc="lower right", frameon=False,
                     fontsize=9, title="percentile", title_fontsize=9)
    leg.get_title().set_color(TEXT_SECONDARY)
    for text in leg.get_texts():
        text.set_color(TEXT_SECONDARY)

    ax.set_xlabel("simulation time t (nondimensional)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_ylabel("shear stress |naive stencil|\n(nondimensional, log scale)",
                   color=TEXT_SECONDARY, fontsize=10)
    ax.tick_params(colors=TEXT_SECONDARY)
    ax.set_title(
        "L10: the 100th percentile separates sharply\nfrom 99th/99.9th, in both codebases",
        color=TEXT_PRIMARY, fontsize=13, fontweight="bold", loc="left", pad=12,
    )
    ax.set_xlim(9.8, 13.6)

    fig.text(
        0.5, 0.005,
        "L10, θ=7°, 32.5 rpm, settled state (past both codebases' own ramp) -- all points from\n"
        "already-completed L10 runs, zero new compute (upstream's OUT_FILES ran for free; ours has one snapshot)",
        ha="center", va="bottom", fontsize=8.5, color=TEXT_SECONDARY, style="italic",
    )

    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(OUT_PATH, dpi=200, facecolor=fig.get_facecolor())
    print(f"\nSaved to {OUT_PATH}")


if __name__ == "__main__":
    main()
