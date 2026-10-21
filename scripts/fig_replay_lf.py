"""Card 04: low fidelity = ONE converged flow period repeated; high fidelity = the normal L8 run (39 conditions, theta 4, 5.5, 7, 9).
Top: ln-log scatter of the repeated-period value against the normal run (line = equality). Bottom: error on unseen conditions against the number of normal runs used.
Data: experiments/multifidelity/replay_lf_l8.csv, replay_lf_mf.json. Usage: uv run python scripts/fig_replay_lf.py"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

QS = (("kLa_1T_10", fs.MK["cstar_10"], r"$k_La$, $C^*\!=\!0.10$ (1/h)"), ("kLa_1T_25", fs.MK["cstar_25"], r"$k_La$, $C^*\!=\!0.25$ (1/h)"), ("dtmix_0.50", fs.MK["chi_0.50"], r"$t_{\mathrm{mix}}$, $\chi\!=\!0.50$ (s)"))


def main():
    rows = list(csv.DictReader((ROOT / "experiments/multifidelity/replay_lf_l8.csv").open()))
    mf = json.loads((ROOT / "experiments/multifidelity/replay_lf_mf.json").read_text())
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.6))
    b = fs.level_colour(8)
    sty = {"MF": (b, "o", "-"), "GPR": ("0.25", "s", "-"), "GPR_LF": ("0.55", "D", "--"), "ratio": (fs.level_colour(9), "^", ":")}
    for j, (q, mk, lab) in enumerate(QS):
        ok = [r for r in rows if r[q] not in ("", "nan") and r[f"hf_{q}"] not in ("", "nan")]
        x = np.array([float(r[f"hf_{q}"]) for r in ok])
        y = np.array([float(r[q]) for r in ok])
        ax = axes[0, j]
        ax.plot(x, y, ls="none", marker=mk, mfc="w", mec=b, ms=6, mew=1.1)
        lim = [min(x.min(), y.min()) * 0.9, max(x.max(), y.max()) * 1.1]
        ax.plot(lim, lim, color=fs.GUIDE, lw=1.0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("normal run: " + lab)
        if j == 0:
            ax.set_ylabel("repeated period")
        ax.grid(**fs.GRID_KW)
        ax = axes[1, j]
        cv = mf[q]["curve"]
        ns = sorted(int(n) for n in cv)
        for k, (col, m, ls) in sty.items():
            ax.plot(ns, [cv[str(n)][k] for n in ns], color=col, marker=m, mfc="w", mec=col, ls=ls, lw=1.0, ms=5)
        ax.set_ylim(0, None)
        ax.set_xlabel("number of normal runs used")
        if j == 0:
            ax.set_ylabel("error on unseen conditions (RMS of ln ratio)")
        ax.grid(**fs.GRID_KW)
    h = [Line2D([], [], color=sty["MF"][0], marker="o", mfc="w", label="mfbml"), Line2D([], [], color=sty["GPR"][0], marker="s", mfc="w", label="GPR, linear mean prior"),
         Line2D([], [], color=sty["GPR_LF"][0], ls="--", marker="D", mfc="w", label="GPR, repeated period as prior"), Line2D([], [], color=sty["ratio"][0], ls=":", marker="^", mfc="w", label="repeated period times a constant")]
    lg = axes[1, 2].legend(handles=h, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / f"experiments/multifidelity/fig_replay_lf.{ext}", dpi=200, bbox_inches="tight", bbox_extra_artists=[lg])


if __name__ == "__main__":
    main()
