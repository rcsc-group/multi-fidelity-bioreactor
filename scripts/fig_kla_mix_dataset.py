"""Figure: every kLa and mixing-time value we have at theta = 7 deg, versus rocking speed, per level.

Panels: kLa_1T_10, kLa_1T_25 (top); dtmix_0.50, dtmix_0.95 (bottom). L8 = all runs (design + 0.625-rpm grid + pairs),
L9 = 10 runs, L10 = the few runs that carry tracer + oxygen (32.5, 35, 37.5 rpm; no L10 rpm sweep exists).
Data: experiments/multifidelity/kla_mix_dataset.csv (scripts/build_kla_mix_dataset.py) + L10 runs listed below.
Usage: uv run python scripts/fig_kla_mix_dataset.py
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

MF = ROOT / "experiments/multifidelity"
L10_RUNS = ("l10c_rpm32.5_seg1", "v_l10_rpm32.5_th7", "fig11_l10_rpm35", "fig11_l10_rpm37.5")
PANELS = (("kLa_1T_10", fs.MK["cstar_10"], r"$k_La$, $C^*\!=\!0.10$ (1/h)"),
          ("kLa_1T_25", fs.MK["cstar_25"], r"$k_La$, $C^*\!=\!0.25$ (1/h)"),
          ("dtmix_0.50", fs.MK["chi_0.50"], r"$t_{\mathrm{mix}}$, $\chi\!=\!0.50$ (s)"),
          ("dtmix_0.95", fs.MK["chi_0.95"], r"$t_{\mathrm{mix}}$, $\chi\!=\!0.95$ (s)"))


def num(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return math.nan
    return v if v > 0 else math.nan


def main():
    rows = [r for r in csv.DictReader((MF / "kla_mix_dataset.csv").open()) if float(r["theta"]) == 7.0]
    pts = {8: [], 9: [], 10: []}
    for r in rows:
        if int(r["level"]) in (8, 9):
            pts[int(r["level"])].append((float(r["rpm"]), {q: num(r[q]) for q, _, _ in PANELS}))
    for run in L10_RUNS:
        d = json.loads((ROOT / "runs" / run / "results.json").read_text())
        om = json.loads((ROOT / "runs" / run / "params.json").read_text())["omega_b"]
        pts[10].append((om * 60 / (2 * math.pi), {q: num(d.get(q)) for q, _, _ in PANELS}))
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.4), sharex=True)
    for ax, (q, mk, lab) in zip(axes.ravel(), PANELS):
        for L in (8, 9, 10):
            xy = np.array([(x, v[q]) for x, v in pts[L] if np.isfinite(v[q])])
            if len(xy):
                kw = fs.series_kw(fs.level_colour(L), mk, stat="mean", ours=True)
                ax.plot(xy[:, 0], xy[:, 1], ls="none", **{k: v for k, v in kw.items() if k != "ls"})
        ax.set_yscale("log")
        ax.set_ylabel(lab)
        ax.grid(**fs.GRID_KW)
    for ax in axes[1]:
        ax.set_xlabel("rocking speed (rpm)")
    handles = [Line2D([], [], ls="none", marker="o", mfc="w", mec=fs.level_colour(L), color=fs.level_colour(L),
                      ms=5, label=f"L{L}") for L in (8, 9, 10)]
    lg = axes[0, 1].legend(handles=handles, frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = MF / f"fig_kla_mix_dataset.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight", bbox_extra_artists=[lg])
        print("wrote", out)
    print({L: len(v) for L, v in pts.items()})


if __name__ == "__main__":
    main()
