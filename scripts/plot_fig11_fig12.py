"""Replicas of Kim et al. (2024) Figs. 11 and 12 -- kLa vs rpm and vs angle.

Both figures come from runs that already exist for Figs 9 and 10: Kim's own
numbers for 9 and 11 are one block of dataset_kimetal2024.xlsx ("Fig. 9,11")
and 10 and 12 are another ("Fig. 10,12"), the same runs scored two ways.
postprocess emits kLa alongside dtmix from every run, so neither figure here
costs a single new cycle.

Estimators map one-to-one onto Kim's columns, and the names were taken from
his: our ``kLa_10/25/50`` is his ``kLa_exp5pts_*`` (a 5-point log-linear fit
of ln(1-C*) around the crossing) and our ``kLa_inst_*`` is his ``kLa_inst_*``
(a finite-difference dC*/dt at the crossing, divided by 1-C*). The panels draw
the 5-point fit, which is the smoother of the two; the printed table carries
both, because where they disagree the oxygen curve is noisy and neither number
should be read as settled.

Marker shape is the saturation threshold and colour is the dataset, the same
encoding as the Fig 9 replica. Levels are never blended: a coarser grid is a
different numerical experiment, not an error bar on this one.

Usage:  uv run python scripts/plot_fig11_fig12.py
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
RUNS = ROOT / "runs"
OUT_DIR = ROOT / "experiments/kimetal2024/figure_replicas"
sys.path.insert(0, str(ROOT))
from scripts.autoextend import tip_results  # noqa: E402
from scripts import figstyle as fs  # noqa: E402

plt.rcParams.update({
    "mathtext.fontset": "cm", "font.family": "serif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
})

C_KIM = fs.KIM
LEVEL_COLOUR = fs.LEVEL_COLOUR
# (our key, Kim's column, marker, label). The C* family is ordered to match
# Fig 9's chi family, so the loosest criterion is "v" in both.
THRESHOLDS = [("kLa_10", "kLa_exp5pts_10", fs.threshold_marker("cstar", 10),
               r"$C^*=10\%$"),
              ("kLa_25", "kLa_exp5pts_25", fs.threshold_marker("cstar", 25),
               r"$C^*=25\%$"),
              ("kLa_50", "kLa_exp5pts_50", fs.threshold_marker("cstar", 50),
               r"$C^*=50\%$")]

FIG11 = {
    "csv": "mixing_kla_vs_frequency.csv", "x": "RPM",
    "values": [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5],
    "series": [(6, "fig9_validate_l6_rpm{v:g}"), (7, "fig9_ladder_l7_rpm{v:g}"),
               (8, "fig9_ladder_l8_rpm{v:g}"), (9, "fig9_l9_rpm{v:g}")],
    "xlabel": r"Rocking frequency $f_b$ (rpm)",
    "title": r"$\theta_{b,max}=7^\circ$", "out": "replicated_Fig11.png",
}
FIG12 = {
    "csv": "mixing_kla_vs_angle.csv", "x": "Angle_deg",
    "values": [2, 3, 4, 5, 6, 7],
    "series": [(9, "fig10_l9_th{v:g}")],
    "xlabel": r"Rocking angle $\theta_{b,max}$ (deg)",
    "title": r"$f_b=32.5$ rpm", "out": "replicated_Fig12.png",
}


def collect(spec: dict) -> pd.DataFrame:
    rows = []
    for level, tmpl in spec["series"]:
        for v in spec["values"]:
            # The chain tip, not the base run: a point finished by a
            # continuation segment keeps NaN in its base results.json.
            d = tip_results(RUNS, tmpl.format(v=v))
            if not d:
                continue
            rows.append({"level": level, "x": float(v),
                         "run_id": tmpl.format(v=v),
                         **{k: d.get(k, math.nan)
                            for k, _, _, _ in THRESHOLDS},
                         **{f"inst_{k}": d.get(f"kLa_inst_{k.split('_')[1]}",
                                               math.nan)
                            for k, _, _, _ in THRESHOLDS}})
    return pd.DataFrame(rows)


def draw(spec: dict) -> None:
    kim = pd.read_csv(ROOT / "experiments/kimetal2024/csv_raw" / spec["csv"])
    kim = kim.sort_values(spec["x"])
    ours = collect(spec)

    print(f"\n{spec['out']}")
    print(f"{'run':<28} {'lvl':>3} {'x':>6} "
          f"{'kLa25':>8} {'inst25':>8} {'Kim25':>8} {'ratio':>7}")
    kim_idx = kim.set_index(spec["x"])
    for _, r in ours.iterrows():
        k25 = (float(kim_idx.loc[r["x"], "kLa_exp5pts_25"])
               if r["x"] in kim_idx.index
               and "kLa_exp5pts_25" in kim_idx.columns else math.nan)
        print(f"{r['run_id']:<28} {r['level']:>3.0f} {r['x']:>6.1f} "
              f"{r['kLa_25']:>8.2f} {r['inst_kLa_25']:>8.2f} "
              f"{k25:>8.2f} {r['kLa_25'] / k25:>7.3f}")

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for _, kim_col, marker, label in THRESHOLDS:
        # Kim's angle block carries only the 10% and 25% crossings; the rpm
        # block carries all three. A missing column is his dataset, not ours.
        if kim_col not in kim.columns:
            continue
        d = kim.dropna(subset=[kim_col])
        ax.plot(d[spec["x"]], d[kim_col],
                **fs.series_kw(C_KIM, marker, stat="max"))
    for level in sorted(ours["level"].unique()) if not ours.empty else []:
        sub = ours[ours["level"] == level]
        colour = LEVEL_COLOUR[int(level)]
        for our_col, _, marker, label in THRESHOLDS:
            d = sub.dropna(subset=[our_col])
            if d.empty:
                continue
            ax.plot(d["x"], d[our_col],
                    **fs.series_kw(colour, marker, stat="max", ours=True))

    # Three small decoder blocks rather than one of 15 entries, and the same
    # three blocks in the same order as Figs 9 and 13: colour names the
    # dataset, marker the threshold, linestyle the provenance.
    from matplotlib.lines import Line2D
    datasets = [Line2D([], [], color=C_KIM, lw=fs.LW_KIM,
                       label="Kim et al.")]
    datasets += [Line2D([], [], color=LEVEL_COLOUR[int(lv)], lw=fs.LW_OURS,
                        ls=":",
                        label=f"L{int(lv)}")
                 for lv in (sorted(ours["level"].unique())
                            if not ours.empty else [])]
    marks = [Line2D([], [], color="0.35", marker=m, ls="none", ms=6, label=lab)
             for _, kim_col, m, lab in THRESHOLDS if kim_col in kim.columns]
    prov = [Line2D([], [], color="0.3", ls="-", lw=fs.LW_KIM,
                   label="Kim et al."),
            Line2D([], [], color="0.3", ls=":", lw=fs.LW_OURS,
                   label="this work")]
    blocks = []
    for handles, title, y in ((datasets, "dataset", 1.0),
                              (marks, "marker", 0.58),
                              (prov, "line", 0.22)):
        lg = ax.legend(handles=handles, fontsize=8, frameon=False,
                       loc="upper left", bbox_to_anchor=(1.02, y),
                       title=title)
        lg.get_title().set_fontsize(8)
        blocks.append(lg)
    for lg in blocks[:-1]:
        ax.add_artist(lg)

    # Log: a coarse grid over-predicts kLa by an order of magnitude (L6 is
    # 10x Kim at 32.5 rpm), and on a linear axis that single point flattens
    # every other series onto the baseline.
    ax.set_yscale("log")
    ax.set_xlabel(spec["xlabel"], fontsize=12)
    ax.set_ylabel(r"$k_La$ (h$^{-1}$)", fontsize=12)
    ax.set_xticks(spec["values"])
    ax.tick_params(which="both", direction="in")
    ax.grid(True, **fs.GRID_KW)
    ax.set_title(spec["title"], loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT_DIR / spec["out"], dpi=150, bbox_inches="tight")
    print(f"saved {OUT_DIR / spec['out']}")


def main() -> None:
    for spec in (FIG11, FIG12):
        draw(spec)


if __name__ == "__main__":
    main()
