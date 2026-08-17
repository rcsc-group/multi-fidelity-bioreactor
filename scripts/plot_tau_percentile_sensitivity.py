"""Figure: shear-stress percentile sensitivity near the top of the tail.

Shows that tau_100 (domain max) is dramatically more sensitive to which
percentile you pick than tau_95/tau_98 are to each other -- the 98th->100th
jump is 9-33x while 95th->98th is only 2.3-3.7x, and growing over time.
Direct evidence for the "one degenerate cut cell dominates the max"
hypothesis (diary.md 2026-08-17): a genuinely heavy-but-smooth tail would
grow steadily across all three percentiles, not blow up only in the last 2%.

Usage:
    uv run python scripts/plot_tau_percentile_sensitivity.py
"""
import numpy as np
import matplotlib.pyplot as plt

SHEAR_STRESS_DAT = "/oscar/scratch/eaguerov/tmp/fork_l10_coldstart/rundir/shear_stress.dat"
OUT_PATH = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/docs_site/assets/img/tau-percentile-sensitivity.png"

# Ordinal blue ramp (dataviz skill palette.md), light->dark for 95th->100th.
COLOR_95 = "#86b6ef"   # step 250
COLOR_98 = "#2a78d6"   # step 450
COLOR_100 = "#104281"  # step 650
COLOR_RATIO_LOW = "#6da7ec"   # step 300 -- 98/95 ratio
COLOR_RATIO_HIGH = "#184f95"  # step 600 -- 100/98 ratio
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e3e2dd"


def load_last_n_rows(path, n=8):
    data = np.loadtxt(path, skiprows=1)
    return data[-n:]


def main():
    data = load_last_n_rows(SHEAR_STRESS_DAT, n=8)
    t = data[:, 1]
    tau_95 = data[:, 2]
    tau_98 = data[:, 3]
    tau_100 = data[:, 4]
    ratio_98_95 = tau_98 / tau_95
    ratio_100_98 = tau_100 / tau_98

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(7.5, 7.2), sharex=True,
        gridspec_kw={"height_ratios": [1.3, 1], "hspace": 0.12},
    )
    fig.patch.set_facecolor("#fcfcfb")

    # --- Panel 1: raw percentiles, log scale ---
    ax1.set_facecolor("#fcfcfb")
    for spine in ("top", "right"):
        ax1.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax1.spines[spine].set_color(GRID)
    ax1.set_yscale("log")
    ax1.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax1.plot(t, tau_100, color=COLOR_100, linewidth=2, marker="o", markersize=5,
              zorder=3, label="100th (max)")
    ax1.plot(t, tau_98, color=COLOR_98, linewidth=2, marker="o", markersize=5,
              zorder=3, label="98th")
    ax1.plot(t, tau_95, color=COLOR_95, linewidth=2, marker="o", markersize=5,
              zorder=3, label="95th")
    # Direct end-of-line labels.
    ax1.annotate("100th (max)", (t[-1], tau_100[-1]), xytext=(6, 0),
                 textcoords="offset points", va="center", color=COLOR_100,
                 fontsize=10, fontweight="bold")
    ax1.annotate("98th", (t[-1], tau_98[-1]), xytext=(6, 0),
                 textcoords="offset points", va="center", color=COLOR_98,
                 fontsize=10, fontweight="bold")
    ax1.annotate("95th", (t[-1], tau_95[-1]), xytext=(6, -2),
                 textcoords="offset points", va="center", color=COLOR_95,
                 fontsize=10, fontweight="bold")
    ax1.set_ylabel("shear stress percentile\n(nondimensional, log scale)",
                    color=TEXT_SECONDARY, fontsize=10)
    ax1.tick_params(colors=TEXT_SECONDARY)
    ax1.set_title(
        "The 100th percentile behaves nothing like the 95th/98th",
        color=TEXT_PRIMARY, fontsize=13, fontweight="bold", loc="left", pad=10,
    )

    # --- Panel 2: ratios between consecutive percentiles, linear scale ---
    ax2.set_facecolor("#fcfcfb")
    for spine in ("top", "right"):
        ax2.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax2.spines[spine].set_color(GRID)
    ax2.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax2.plot(t, ratio_100_98, color=COLOR_RATIO_HIGH, linewidth=2, marker="o",
              markersize=5, zorder=3)
    ax2.plot(t, ratio_98_95, color=COLOR_RATIO_LOW, linewidth=2, marker="o",
              markersize=5, zorder=3)
    ax2.annotate("100th / 98th ratio", (t[-1], ratio_100_98[-1]), xytext=(6, 0),
                 textcoords="offset points", va="center", color=COLOR_RATIO_HIGH,
                 fontsize=10, fontweight="bold")
    ax2.annotate("98th / 95th ratio", (t[-1], ratio_98_95[-1]), xytext=(6, 0),
                 textcoords="offset points", va="center", color=COLOR_RATIO_LOW,
                 fontsize=10, fontweight="bold")
    ax2.set_ylabel("ratio between\nconsecutive percentiles", color=TEXT_SECONDARY,
                    fontsize=10)
    ax2.set_xlabel("simulation time t (nondimensional)", color=TEXT_SECONDARY,
                    fontsize=10)
    ax2.tick_params(colors=TEXT_SECONDARY)
    ax2.set_ylim(0, max(ratio_100_98) * 1.25)

    fig.text(
        0.5, 0.005,
        "L10, θ=7°, 32.5 rpm, settled state -- 100th/98th jumps 9-33x (and growing) "
        "while 98th/95th stays flat at 2.3-3.7x",
        ha="center", va="bottom", fontsize=9, color=TEXT_SECONDARY, style="italic",
    )

    for ax in (ax1, ax2):
        ax.set_xlim(t[0] - 0.005, t[-1] + 0.03)

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT_PATH, dpi=200, facecolor=fig.get_facecolor())
    print(f"Saved to {OUT_PATH}")

    print("\n=== Table version ===")
    print(f"{'t':>6} {'tau_95':>10} {'tau_98':>10} {'tau_100':>10} {'100/98':>8} {'98/95':>8}")
    for i in range(len(t)):
        print(f"{t[i]:>6.2f} {tau_95[i]:>10.6g} {tau_98[i]:>10.6g} {tau_100[i]:>10.6g} "
              f"{ratio_100_98[i]:>8.1f} {ratio_98_95[i]:>8.2f}")


if __name__ == "__main__":
    main()
