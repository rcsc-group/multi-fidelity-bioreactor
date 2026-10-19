"""Figure: L8 -> L9 prediction of kLa_1T_10 and dtmix_0.50 at theta 7 (test S89, experiments/multifidelity/test_s89_l8_to_l9.json).
Blue dashed: L8 fit. Green: L9 runs. Grey bars: prediction (95% range) for each L9 run made WITHOUT any L9 run within 0.25 rpm of it.
Usage: uv run python scripts/fig_s89.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

LAB = {"kLa_1T_10": (fs.MK["cstar_10"], r"$k_La$, $C^*\!=\!0.10$ (1/h)"), "dtmix_0.50": (fs.MK["chi_0.50"], r"$t_{\mathrm{mix}}$, $\chi\!=\!0.50$ (s)")}


def main():
    d = json.loads((ROOT / "experiments/multifidelity/test_s89_l8_to_l9.json").read_text())
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.9))
    for ax, q in zip(axes, LAB):
        rows = sorted(d[q]["rows"], key=lambda r: r["rpm"])
        mk = LAB[q][0]
        x = [r["rpm"] for r in rows]
        ax.plot(x, [math.exp(r["l8"]) for r in rows], color=fs.level_colour(8), ls="--", lw=1.0)
        for r in rows:
            m, lo, hi = math.exp(r["mf"]), math.exp(r["lo"]), math.exp(r["hi"])
            ax.errorbar(r["rpm"] + 0.3, m, yerr=[[m - lo], [hi - m]], color="0.25", marker="_", ms=6, elinewidth=1.0, capsize=2, ls="none")
        kw = fs.series_kw(fs.level_colour(9), mk, stat="mean", ours=True)
        kw.pop("ls")
        ax.plot(x, [math.exp(r["obs"]) for r in rows], ls="none", **kw)
        ax.set_yscale("log")
        ax.set_xlabel("rocking speed (rpm)")
        ax.set_ylabel(LAB[q][1])
        ax.grid(**fs.GRID_KW)
    h = [Line2D([], [], color=fs.level_colour(8), ls="--", label="L8 fit"),
         Line2D([], [], ls="none", marker="o", mfc="w", mec=fs.level_colour(9), label="L9"),
         Line2D([], [], color="0.25", marker="_", ls="none", ms=8, label="prediction without nearby L9 runs")]
    lg = axes[1].legend(handles=h, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / f"experiments/multifidelity/fig_s89.{ext}", dpi=200, bbox_inches="tight", bbox_extra_artists=[lg])


if __name__ == "__main__":
    main()
