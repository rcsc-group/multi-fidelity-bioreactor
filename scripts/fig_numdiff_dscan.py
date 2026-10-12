"""Figure: kLa and mixing time do not respond to the physical scalar diffusivity at L7-L9 (test K), and their value
is set by the grid level (L7 -> L10 at physical D). One figure, two panels, our own solver only.

Data: experiments/multifidelity/test_k_summary.json (k_l<L>_s<scale> runs, 32.5 rpm, theta 7, D scaled by s for
tracer and oxygen in the liquid) and runs/l10c_rpm32.5_seg1/results.json (cold L10, physical D, s = 1).
Reference line: kLa proportional to D^(1/2) (penetration theory, resolved interfacial layer), anchored at the L9 s = 1
value. Usage: uv run python scripts/fig_numdiff_dscan.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

LEVELS = (7, 8, 9)
PANELS = (("kLa_1T_25", r"$k_La$ (1/h)", fs.MK["cstar_25"]),
          ("dtmix_0.95", r"$\Delta t_{\mathrm{mix},0.95}$ (s)", fs.MK["chi_0.95"]))


def main():
    runs = json.loads((ROOT / "experiments/multifidelity/test_k_summary.json").read_text())["runs"]
    l10 = json.loads((ROOT / "runs/l10c_rpm32.5_seg1/results.json").read_text())
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), sharex=True)
    for ax, (q, ylab, mk) in zip(axes, PANELS):
        for L in LEVELS:
            pts = sorted((float(s), v[q]) for s, v in runs[f"L{L}"].items() if v.get(q) is not None)
            s, y = np.array(pts).T
            ax.plot(s, y, label=f"L{L}", **fs.series_kw(fs.level_colour(L), mk, stat="mean", ours=True))
        ax.plot([1.0], [l10[q]], label="L10", **fs.series_kw(fs.level_colour(10), mk, stat="mean", ours=True))
        if q == "kLa_1T_25":
            y9 = runs["L9"]["1e0"][q]
            sr = np.logspace(0, 3, 20)
            ax.plot(sr, y9 * sr ** 0.5, ls="--", color=fs.NEUTRAL, lw=1.0, label=r"$\propto D^{1/2}$")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$D/D_{\mathrm{phys}}$")
        ax.set_ylabel(ylab)
        ax.grid(**fs.GRID_KW)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = ROOT / f"experiments/multifidelity/fig_numdiff_dscan.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print("wrote", out)


if __name__ == "__main__":
    main()
