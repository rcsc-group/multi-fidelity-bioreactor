"""Figure: H3 (L8 -> L9 -> L10 KRR-LR-GPR on log QoI) leave-one-out prediction vs the L10 simulation, 13 points.

(a) tau_95, (b) EDR. Error bars: MF 95% interval. Grey: the GP trained on the L10 points only (same LOO).
Solid line: prediction = simulation; dashed: +-10% (the pre-registered RMSE bar).
Data: experiments/multifidelity/test_h3_mfbml.json. Usage: uv run python scripts/fig_h3_parity.py
"""
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

MF = ROOT / "experiments/multifidelity"
PANELS = (("tau_95_t", fs.MK["tau"], "max", r"$\tau_{95}$ (Pa)"),
          ("edr_mean_t", fs.MK["ediss"], "mean", r"EDR (W/m$^3$)"))


def main():
    res = json.loads((MF / "test_h3_mfbml.json").read_text())["results"]
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8))
    c10 = fs.level_colour(10)
    for ax, (q, mk, stat, lab) in zip(axes, PANELS):
        p = res[q]["points"]
        obs = np.array([x["obs"] for x in p])
        mf = np.array([x["mf_pred"] for x in p])
        err = np.array([[x["mf_pred"] - x["mf_lo"] for x in p], [x["mf_hi"] - x["mf_pred"] for x in p]])
        hf = np.array([x["hf_pred"] for x in p])
        ax.plot(obs, hf, ls="none", **{k: v for k, v in fs.series_kw(fs.NEUTRAL, mk, stat=stat).items() if k != "ls"})
        ax.errorbar(obs, mf, yerr=err, ls="none", ecolor=c10, elinewidth=0.8, capsize=2,
                    **{k: v for k, v in fs.series_kw(c10, mk, stat=stat, ours=True).items() if k != "ls"})
        lo, hi = obs.min() / 1.6, obs.max() * 1.6
        xs = np.geomspace(lo, hi, 20)
        ax.plot(xs, xs, color=fs.NEUTRAL, lw=1.0)
        for f in (0.9, 1.1):
            ax.plot(xs, f * xs, color=fs.NEUTRAL, lw=0.8, ls="--")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_xlabel(f"{lab}, L10 simulation")
        ax.set_ylabel(f"{lab}, prediction")
        ax.grid(**fs.GRID_KW)
        s = res[q]["stats"]
        print(f"{q}: MF {100 * s['mf_rmse']:.1f}%  L10-only GP {100 * s['hf_rmse']:.1f}%")
    handles = [Line2D([], [], color=c10, lw=4, label="L8 → L9 → L10 (mfbml)"),
               Line2D([], [], color=fs.NEUTRAL, lw=4, label="L10 only (GP)"),
               Line2D([], [], color=fs.NEUTRAL, lw=1.0, label="exact"),
               Line2D([], [], color=fs.NEUTRAL, lw=0.8, ls="--", label=r"$\pm$10%")]
    lg = axes[1].legend(handles=handles, fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = MF / f"fig_h3_parity.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight", bbox_extra_artists=[lg])
        print("wrote", out)


if __name__ == "__main__":
    main()
