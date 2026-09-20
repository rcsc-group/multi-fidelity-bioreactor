"""Replica of Kim et al. (2024) Fig. 9(a) -- mixing time vs rocking frequency.

Kim's panel (a): mixing times at chi = 0.95 (circles), 0.75 (upper
triangles) and 0.50 (lower triangles) on the left axis, and the spatially
averaged absolute steady-streaming vorticity in the water (hollow squares)
on the right axis, for theta_max = 7 deg. Panels (b)-(d) are tracer-field
snapshots at chi = 0.5 for three rpm; not reproduced here.

Reference values are Kim's OWN published numbers, not a digitization:
csv_raw/mixing_kla_vs_frequency.csv is a bit-exact export of
dataset_kimetal2024.xlsx sheet "Fig. 9,11" (see csv_raw/PROVENANCE.md).

Our values come from results.json, computed by postprocess's chi -> dtmix
pipeline, which returns NaN unless sigma^2_max is within 20% of its analytic
0.25 (tests/verification/test_mixing_metrics.py). A level whose points are
missing therefore failed that invariant rather than simply not having run --
the printed table says which.

KIM'S RESULTS ARE AT n_L = 2^10, i.e. our L10 (Main.tex sec. 4: "A uniform
mesh with n_L=2^10 grid cells along the width of the domain is used. Each
grid cell has a dimensional size of 0.24 mm, resulting in a total of
1.05e6 cells"). Our levels map onto his one-to-one -- our domain width is
L_bio = L_x = 0.25 m with 2^L cells across it, so L10 gives 0.244 mm and
1024^2 = 1.05e6 cells, matching both of his stated numbers.

So a coarser level sitting below Kim's curve is NOT a discrepancy; it is a
grid 2^(10-L) times coarser per direction, and tracer numerical diffusion is
the quantity most sensitive to cell size. Measured: L6 12.9x fast, L7 6.5x
on dtmix_0.95, a per-level gain of 1.98 that extrapolates to parity at
~L10 -- Kim's own mesh (diary.md 2026-09-15 (9)).

Levels are therefore plotted separately and never blended, and the L-vs-Kim
ratio is only a convergence statement, never an error.

Usage:  uv run python scripts/plot_fig9.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
KIM_CSV = ROOT / "experiments/kimetal2024/csv_raw/mixing_kla_vs_frequency.csv"
OUT = ROOT / "experiments/kimetal2024/figure_replicas/replicated_Fig9.png"
RUNS = ROOT / "runs"
sys.path.insert(0, str(ROOT))
from scripts.autoextend import tip_results  # noqa: E402
from scripts import figstyle as fs  # noqa: E402

RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
# (level, run_id template, colour). Ladder points are single-rpm probes;
# sweep points use the fig9_ prefix.
SERIES = [
    (6, "fig9_validate_l6_rpm{rpm:g}", fs.level_colour(6)),
    (7, "fig9_l7_rpm{rpm:g}", fs.level_colour(7)),
    (7, "fig9_ladder_l7_rpm{rpm:g}", fs.level_colour(7)),
    (8, "fig9_ladder_l8_rpm{rpm:g}", fs.level_colour(8)),
    (9, "fig9_ladder_l9_rpm{rpm:g}", fs.level_colour(9)),
    (9, "fig9_l9_rpm{rpm:g}", fs.level_colour(9)),   # the sweep itself
]
# Marker ranks the homogeneity criterion (v -> ^ -> P with stringency) and is
# shared with Fig 10; it never collides with a physical-quantity shape.
THRESHOLDS = [("dtmix_0.95", fs.threshold_marker("chi", 0.95), 0.95),
              ("dtmix_0.75", fs.threshold_marker("chi", 0.75), 0.75),
              ("dtmix_0.50", fs.threshold_marker("chi", 0.50), 0.50)]


def collect() -> pd.DataFrame:
    rows = []
    for level, tmpl, colour in SERIES:
        for rpm in RPMS:
            # A point that ran out of clock is finished by a continuation
            # segment (scripts/autoextend.py); its completed measurement lives
            # in the LAST segment, while the base run's own results.json keeps
            # its NaN forever. Reading the base run_id directly would report
            # every extended point as missing.
            d = tip_results(RUNS, tmpl.format(rpm=rpm))
            if not d:
                continue
            rows.append({
                "level": level, "rpm": float(rpm), "colour": colour,
                "run_id": tmpl.format(rpm=rpm),
                "sigma2_max": d.get("sigma2_max", math.nan),
                # Kim's right-hand axis is the STEADY-STREAMING vorticity (curl of the
                # time-averaged flow), not vor_mean (time-average of |curl u|).
                # See postprocess docstring; vor_mean is NaN-free but wrong here.
                "vor_mean": d.get("vor_streaming", math.nan),
                **{k: d.get(k, math.nan) for k, _, _ in THRESHOLDS},
            })
    return pd.DataFrame(rows).sort_values(["level", "rpm"])


def main() -> None:
    kim = pd.read_csv(KIM_CSV).sort_values("RPM")
    ours = collect()

    print(f"{'run':<28} {'lvl':>3} {'rpm':>6} {'s2max':>7} "
          f"{'dt50':>8} {'dt95':>8} {'dt95/Kim':>9}")
    kim_idx = kim.set_index("RPM")
    for _, r in ours.iterrows():
        ratio = (r["dtmix_0.95"] / kim_idx.loc[r["rpm"], "dtmix_strict_0.95"]
                 if not math.isnan(r["dtmix_0.95"]) else math.nan)
        print(f"{r['run_id']:<28} {r['level']:>3} {r['rpm']:>6.1f} "
              f"{r['sigma2_max']:>7.4f} {r['dtmix_0.50']:>8.2f} "
              f"{r['dtmix_0.95']:>8.2f} {ratio:>9.3f}")
    if ours.empty:
        print("no fig9 runs with results.json yet")
        return

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax2 = ax.twinx()

    for col, marker, thr in THRESHOLDS:
        ax.plot(kim["RPM"], kim[f"dtmix_strict_{thr:g}"],
                **fs.series_kw(fs.KIM, marker, stat="max"),
                label=rf"Kim  $\chi={thr:.2f}$")
    # Steady-streaming vorticity is a time-and-space mean, so it is hollow.
    ax2.plot(kim["RPM"], kim["vor_meanabs_steady_streaming"],
             **fs.series_kw(fs.KIM, fs.MK["vorticity"], stat="mean"),
             label=r"Kim  $\langle|\bar\xi_b'|\rangle$")

    for level in sorted(ours["level"].unique()):
        sub = ours[ours["level"] == level]
        colour = sub["colour"].iloc[0]
        for col, marker, thr in THRESHOLDS:
            d = sub.dropna(subset=[col])
            if d.empty:
                continue
            ax.plot(d["rpm"], d[col],
                    **fs.series_kw(colour, marker, stat="max", ours=True),
                    label=rf"L{level}  $\chi={thr:.2f}$")
        dv = sub.dropna(subset=["vor_mean"])
        if not dv.empty:
            ax2.plot(dv["rpm"], dv["vor_mean"],
                     **fs.series_kw(colour, fs.MK["vorticity"], stat="mean",
                                    ours=True),
                     label=rf"L{level}  $\langle|\bar\xi_b'|\rangle$")

    ax.set_yscale("log")
    ax.set_xlabel(r"Rocking frequency $f_b$ (rpm)", fontsize=12)
    ax.set_ylabel("Mixing time (s)", fontsize=12)
    ax2.set_ylabel(r"$\langle|\bar\xi_b'|\rangle$ (1/s)", fontsize=12)
    ax.tick_params(which="both", direction="in")
    ax2.tick_params(which="both", direction="in")
    ax.grid(True, **fs.GRID_KW)
    ax.set_xticks(RPMS)
    ax.set_title(r"$\theta_{b,max}=7^\circ$", loc="left", fontsize=10)

    # Three small decoder blocks rather than one enumerated list, matching
    # Fig 13. Enumerating every (dataset x threshold) pair costs fourteen
    # entries to say what three keys say once each, and it hides the fact
    # that the encoding is the same in both figures.
    from matplotlib.lines import Line2D
    levels_shown = sorted(ours["level"].unique())
    ds = [Line2D([], [], color=fs.KIM, lw=fs.LW_KIM, label="Kim et al.")]
    ds += [Line2D([], [], color=fs.level_colour(lv), lw=fs.LW_OURS, ls=":",
                  label=f"L{int(lv)}") for lv in levels_shown]
    enc = [Line2D([], [], color="0.3", marker=m, ls="none", mfc="0.3", ms=6,
                  label=rf"$\chi={thr:.2f}$")
           for _, m, thr in THRESHOLDS]
    enc.append(Line2D([], [], color="0.3", marker=fs.MK["vorticity"],
                      ls="none", mfc="w", mec="0.3", ms=6,
                      label=r"$\langle|\bar\xi_b'|\rangle$  (hollow: a mean)"))
    prov = [Line2D([], [], color="0.3", ls="-", lw=fs.LW_KIM,
                   label="Kim et al."),
            Line2D([], [], color="0.3", ls=":", lw=fs.LW_OURS,
                   label="this work")]
    legs = []
    for handles, title, y in ((ds, "dataset", 0.95), (enc, "marker / fill", 0.58),
                              (prov, "line", 0.22)):
        lg = fig.legend(handles=handles, fontsize=8.5, loc="upper left",
                        bbox_to_anchor=(1.0, y), frameon=False, title=title)
        lg.get_title().set_fontsize(8.5)
        legs.append(lg)
    for lg in legs[:-1]:
        fig.add_artist(lg)
    fig.tight_layout()
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"\nsaved {OUT}")


if __name__ == "__main__":
    main()
