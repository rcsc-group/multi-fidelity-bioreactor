"""Replica of Kim et al. (2024) Fig. 10(a) -- mixing time vs rocking angle.

The angle counterpart of Fig. 9: mixing time on the left axis and the
spatially averaged absolute steady-streaming vorticity on the right, against
rocking angle at f_b = 32.5 rpm.

Kim publishes only chi = 0.50 for the angle series -- his Fig. 9 carries 0.95,
0.75 and 0.50, his Fig. 10 carries 0.50 alone -- so that is the only threshold
with a reference to compare against, and it is the one the angle sweep targets
(scripts/submit_fig10.py). The tighter thresholds are still computed and are
printed in the table; they simply have nothing to sit beside.

Reference values are Kim's own published numbers, a bit-exact export of the
"Fig. 10,12" block of his spreadsheet, not a digitisation.

As in Fig. 9, levels are plotted separately and never blended: a coarser grid
is a different numerical experiment, and tracer numerical diffusion is the
quantity most sensitive to cell size, so a level sitting below Kim's curve is
a convergence statement rather than an error.

Usage:  uv run python scripts/plot_fig10.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
KIM_CSV = ROOT / "experiments/kimetal2024/csv_raw/mixing_kla_vs_angle.csv"
OUT = ROOT / "experiments/kimetal2024/figure_replicas/replicated_Fig10.png"
RUNS = ROOT / "runs"
sys.path.insert(0, str(ROOT))
from scripts.autoextend import tip_results  # noqa: E402

plt.rcParams.update({
    "mathtext.fontset": "cm", "font.family": "serif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
})

ANGLES = [2, 3, 4, 5, 6, 7]
SERIES = [(9, "fig10_l9_th{v:g}", "#CC79A7")]
# Kim's angle block has a reference for chi = 0.50 only; the others are ours
# alone and are printed rather than plotted.
PLOTTED = [("dtmix_0.50", "v", 0.50, "dtmix_strict_0.5")]
EXTRA = ["dtmix_0.75", "dtmix_0.95"]


def collect() -> pd.DataFrame:
    rows = []
    for level, tmpl, colour in SERIES:
        for v in ANGLES:
            d = tip_results(RUNS, tmpl.format(v=v))
            if not d:
                continue
            rows.append({"level": level, "theta": float(v), "colour": colour,
                         "run_id": tmpl.format(v=v),
                         "sigma2_max": d.get("sigma2_max", math.nan),
                         "vor": d.get("vor_streaming", math.nan),
                         **{k: d.get(k, math.nan)
                            for k in [p[0] for p in PLOTTED] + EXTRA}})
    # An empty frame has no columns to sort by, and this figure is expected to
    # be empty until the angle sweep lands.
    df = pd.DataFrame(rows)
    return df.sort_values(["level", "theta"]) if not df.empty else df


def main() -> None:
    kim = pd.read_csv(KIM_CSV).sort_values("Angle_deg")
    ours = collect()

    print(f"{'run':<24} {'lvl':>3} {'deg':>5} {'s2max':>7} "
          f"{'dt50':>9} {'dt75':>9} {'dt95':>9} {'dt50/Kim':>9}")
    kim_idx = kim.set_index("Angle_deg")
    for _, r in ours.iterrows():
        ref = float(kim_idx.loc[r["theta"], "dtmix_strict_0.5"])
        ratio = r["dtmix_0.50"] / ref if ref else math.nan
        print(f"{r['run_id']:<24} {r['level']:>3.0f} {r['theta']:>5.1f} "
              f"{r['sigma2_max']:>7.4f} {r['dtmix_0.50']:>9.2f} "
              f"{r['dtmix_0.75']:>9.2f} {r['dtmix_0.95']:>9.2f} {ratio:>9.3f}")
    if ours.empty:
        print("no fig10 runs with results yet -- submit scripts/submit_fig10.py")

    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax2 = ax.twinx()
    for col, marker, thr, kim_col in PLOTTED:
        ax.plot(kim["Angle_deg"], kim[kim_col], color="black", marker=marker,
                ms=6, lw=1.3, label=rf"Kim  $\chi={thr:.2f}$")
    ax2.plot(kim["Angle_deg"], kim["vor_meanabs_steady_streaming"],
             color="#E69F00", marker="s", ms=6, lw=1.3, ls="--", mfc="w",
             label=r"Kim  $\langle|\bar\xi_b'|\rangle$")

    for level in sorted(ours["level"].unique()) if not ours.empty else []:
        sub = ours[ours["level"] == level]
        colour = sub["colour"].iloc[0]
        for col, marker, thr, _ in PLOTTED:
            d = sub.dropna(subset=[col])
            if d.empty:
                continue
            ax.plot(d["theta"], d[col], color=colour, marker=marker, ms=5,
                    lw=1.1, ls=":", label=rf"L{int(level)}  $\chi={thr:.2f}$")
        dv = sub.dropna(subset=["vor"])
        if not dv.empty:
            ax2.plot(dv["theta"], dv["vor"], color=colour, marker="s", ms=4,
                     lw=1.0, ls="-.", mfc="w",
                     label=rf"L{int(level)}  $\langle|\bar\xi_b'|\rangle$")

    ax.set_yscale("log")
    ax.set_xlabel(r"Rocking angle $\theta_{b,max}$ (deg)", fontsize=12)
    ax.set_ylabel("Mixing time (s)", fontsize=12)
    ax2.set_ylabel(r"$\langle|\bar\xi_b'|\rangle$ (1/s)", fontsize=12,
                   color="#E69F00")
    ax2.tick_params(axis="y", colors="#E69F00")
    ax.tick_params(which="both", direction="in")
    ax2.tick_params(which="both", direction="in")
    ax.grid(True, which="major", ls=":", alpha=0.4)
    ax.set_xticks(ANGLES)
    ax.set_title(r"$f_b=32.5$ rpm", loc="left", fontsize=10)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    fig.legend(h1 + h2, l1 + l2, fontsize=8, loc="center left",
               bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"\nsaved {OUT}")


if __name__ == "__main__":
    main()
