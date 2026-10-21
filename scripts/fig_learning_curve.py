"""Figure: error of the L10 prediction against the cost of the training data (core-hours), hydrodynamic quantities.
From experiments/multifidelity/h3_baselines_learning_curve.json. Cost per training run: L10 935 core-h; L9 runs at all 13 design points: 13 x 74 core-h (diary 2026-10-11/13 costs).
L10 only: n L10 runs. L9 plus a ratio, and the three-level chain: n L10 runs + L9 at all 13 points.
Usage: uv run python scripts/fig_learning_curve.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

C10, C9, L9_COST = 935.0, 74.0, 13 * 74.0
LAB = {"tau_95_t": r"$\tau_{95}$", "tau_mean_t": r"$\overline{\tau}$", "edr_mean_t": "EDR"}


def main():
    d = json.loads((ROOT / "experiments/multifidelity/h3_baselines_learning_curve.json").read_text())
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), sharey=True)
    for ax, q in zip(axes, LAB):
        lc = {int(n): v for n, v in d[q]["learning_curve"].items()}
        ns = sorted(lc)
        for key, colour, ls, mk, cost in (("l10_only", fs.level_colour(10), "-", "o", lambda n: n * C10),
                                          ("l9_loglin", fs.level_colour(9), ":", "s", lambda n: n * C10 + L9_COST),
                                          ("chain", "0.25", "-", "^", lambda n: n * C10 + L9_COST)):
            pts = [(cost(n), 100 * lc[n][key]) for n in ns if key in lc[n]]
            ax.plot(*zip(*pts), color=colour, ls=ls, marker=mk, mfc="w", mec=colour, ms=5, lw=1.0)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_title(LAB[q], fontsize=10)
        ax.set_xlabel("training cost (core-hours)")
        ax.set_xticks([3000, 5000, 10000])
        ax.set_xticklabels(["3k", "5k", "10k"])
        ax.minorticks_off()
        ax.grid(**fs.GRID_KW)
    axes[0].set_ylabel("error on unseen L10 runs (%)")
    h = [Line2D([], [], color=fs.level_colour(10), marker="o", mfc="w", label="L10 runs only"),
         Line2D([], [], color=fs.level_colour(9), ls=":", marker="s", mfc="w", label="L9 times a fitted ratio"),
         Line2D([], [], color="0.25", marker="^", mfc="w", label="three-level chain")]
    lg = axes[2].legend(handles=h, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / f"experiments/multifidelity/fig_learning_curve.{ext}", dpi=200, bbox_inches="tight", bbox_extra_artists=[lg])


if __name__ == "__main__":
    main()
