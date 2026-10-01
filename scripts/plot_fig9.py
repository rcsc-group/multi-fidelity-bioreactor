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

Drawing (2026-10-01): Kim's own encoding, not the project's level colours --
colour+marker name the threshold, hollow purple squares the vorticity, log-log
axes over his ranges. Only the L9 sweep is drawn and named in the key; Kim is
markers only, L9 the same markers joined by a dotted line. The printed table
still lists every level.

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
OUT = ROOT / "experiments/kimetal2024/figure_replicas/replicated_Fig09.png"
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
    # L10 kLa points (submit_fig11_l10_kla.py) carry the tracer too. They
    # release 5 cycles after a warm restart, not at absolute cycle 80.
    (10, "fig11_l10_rpm{rpm:g}", fs.level_colour(10)),
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

    # Kim's Fig 9(a) encoding, exactly (Figures/Fig_dtmix_rpm.pdf): colour and
    # marker name the THRESHOLD -- chi=0.95 blue circles, 0.75 red up-triangles,
    # 0.50 green down-triangles -- and the vorticity is hollow purple squares on
    # a purple right axis. Log-log, f_b 10-100 rpm, dtmix 1-1000 s, vorticity
    # 0.1-10. Kim draws markers only; ours repeat his markers joined by a
    # dotted line, which is what tells the two datasets apart.
    # One level only (L9, the full sweep), named in the key.
    sweep = ours[ours["run_id"].str.startswith("fig9_l9_rpm")]
    enc = {0.95: ("#1f5fa8", "o"), 0.75: ("#d62728", "^"), 0.50: ("#2ca02c", "v")}
    purple = "#9b4f9b"
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    ax2 = ax.twinx()
    for col, _, thr in THRESHOLDS:
        c, m = enc[thr]
        ax.plot(kim["RPM"], kim[f"dtmix_strict_{thr:g}"], ls="none", marker=m,
                color=c, ms=7)
        d = sweep.dropna(subset=[col])
        ax.plot(d["rpm"], d[col], ls=":", lw=1.1, marker=m, color=c, ms=5.5)
    ax2.plot(kim["RPM"], kim["vor_meanabs_steady_streaming"], ls="none",
             marker="s", mfc="w", mec=purple, mew=1.3, ms=7)
    dv = sweep.dropna(subset=["vor_mean"])
    ax2.plot(dv["rpm"], dv["vor_mean"], ls=":", lw=1.1, color=purple,
             marker="s", mfc="w", mec=purple, mew=1.1, ms=5.5)

    for a_ in (ax, ax2):
        a_.set_xscale("log")
        a_.set_yscale("log")
        a_.tick_params(which="both", direction="in")
    from matplotlib.ticker import NullFormatter
    ax.xaxis.set_minor_formatter(NullFormatter())   # decades only, as Kim
    ax.set_xlim(10, 100)
    ax.set_ylim(1, 1000)
    ax2.set_ylim(0.1, 10)
    ax.set_xlabel(r"$f_b$ (rpm)", fontsize=12)
    ax.set_ylabel(r"$\Delta t_{\mathrm{mix}}$ (s)", fontsize=12)
    ax2.set_ylabel(r"$\langle|\overline{\xi_b'}|\rangle$ (1/s)", fontsize=12,
                   color=purple)
    ax2.tick_params(axis="y", which="both", colors=purple)
    ax2.spines["right"].set_color(purple)

    from matplotlib.lines import Line2D
    quantity = [Line2D([], [], ls="none", marker=enc[t][1], color=enc[t][0],
                       ms=6.5, label=rf"$\chi={t:.2f}$") for _, _, t in THRESHOLDS]
    quantity.append(Line2D([], [], ls="none", marker="s", mfc="w", mec=purple,
                           mew=1.3, ms=6.5,
                           label=r"$\langle|\overline{\xi_b'}|\rangle$"))
    dataset = [Line2D([], [], ls="none", marker="o", color="0.3", ms=6.5,
                      label="Kim et al."),
               Line2D([], [], ls=":", lw=1.1, marker="o", color="0.3", ms=5,
                      label="L9")]
    lg1 = fig.legend(handles=quantity, fontsize=9, frameon=False,
                     loc="upper left", bbox_to_anchor=(1.0, 0.95))
    fig.add_artist(lg1)
    fig.legend(handles=dataset, fontsize=9, frameon=False,
               loc="upper left", bbox_to_anchor=(1.0, 0.55))
    fig.tight_layout()
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"\nsaved {OUT}")


if __name__ == "__main__":
    main()
