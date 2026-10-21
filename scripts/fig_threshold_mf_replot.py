"""Redraw fig_threshold_mf from experiments/multifidelity/threshold_mf.json (no recomputation). Usage: uv run python scripts/fig_threshold_mf_replot.py"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402
from scripts import postprocess as pp  # noqa: E402
from scripts.fig_threshold_mf import PAIRS, SIZES  # noqa: E402


def main():
    out = json.loads((ROOT / "experiments/multifidelity/threshold_mf.json").read_text())
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(1, 5, figsize=(17.0, 3.5))
    ax = axes[0]
    run = ROOT / "runs" / "kmix_l8_rpm25"
    p = json.loads((run / "params.json").read_text())
    t, c = pp._compute_c_star(run)
    T_bio, T_nd = pp._t_scales(p)
    tt = (t - t[np.argmax(c > 1e-6)]) * T_bio
    ax.plot(tt, c, color=fs.level_colour(8), lw=1.0)
    for x in (0.10, 0.25, 0.50, 0.75):
        ax.axhline(x, color=fs.GUIDE, lw=1.0, ls=":")
    ax.set_xlim(0, 360)
    ax.set_xlabel("time after release (s)")
    ax.set_ylabel(r"$C^*$")
    ax.set_title("the problem: L8, 25 rpm, 7$^\\circ$", fontsize=10)
    ax.grid(**fs.GRID_KW)
    sty = {"MF": (fs.level_colour(8), "o", "-"), "GPR": ("0.25", "s", "-"), "GPR_LF": ("0.55", "D", "--"), "ratio": (fs.level_colour(9), "^", ":")}
    for ax, (lo, hi) in zip(axes[1:], PAIRS):
        cv = out[f"{lo}->{hi}"]
        for k, (col, mk, ls) in sty.items():
            ax.plot(SIZES, [cv[str(n)][k] for n in SIZES], color=col, marker=mk, mfc="w", mec=col, ls=ls, lw=1.0, ms=5)
        ax.set_xscale("log")
        ax.minorticks_off()
        ax.set_ylim(0, None)
        ax.set_xticks([4, 8, 16, 32])
        ax.set_xticklabels(["4", "8", "16", "32"])
        ax.set_xlabel("runs with the higher threshold")
        ax.set_title(rf"$C^*$ {lo / 100:.2f} $\to$ {hi / 100:.2f}", fontsize=10)
        ax.grid(**fs.GRID_KW)
    axes[1].set_ylabel("error on unseen runs (RMS of ln ratio)")
    h = [Line2D([], [], color=sty["MF"][0], marker="o", mfc="w", label="mfbml"), Line2D([], [], color=sty["GPR"][0], marker="s", mfc="w", label="GPR, linear mean prior"),
         Line2D([], [], color=sty["GPR_LF"][0], ls="--", marker="D", mfc="w", label="GPR, lower threshold as prior"),
         Line2D([], [], color=sty["ratio"][0], ls=":", marker="^", mfc="w", label="lower threshold times a constant")]
    lg = axes[-1].legend(handles=h, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / f"experiments/multifidelity/fig_threshold_mf.{ext}", dpi=200, bbox_inches="tight", bbox_extra_artists=[lg])


if __name__ == "__main__":
    main()
