"""Figure: how well H3 (L8 -> L9 -> L10 KRR-LR-GPR on log QoI) predicts each L10 QoI, leave-one-out, 13 points.

One panel per QoI: tau_95, tau_mean, EDR. Bars: 95% predictive interval. Solid: prediction = L10 simulation;
dashed: +-10% (pre-registered RMSE bar). Panel title: the QoI and its LOO relative RMSE.
Data: experiments/multifidelity/test_h3_mfbml.json. Usage: uv run python scripts/fig_h3_parity.py
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

MF = ROOT / "experiments/multifidelity"
PANELS = (("tau_95_t", fs.MK["tau"], "max", r"$\tau_{95}$", "Pa"),
          ("tau_mean_t", fs.MK["tau"], "mean", r"$\overline{\tau}$", "Pa"),
          ("edr_mean_t", fs.MK["ediss"], "mean", "EDR", r"W/m$^3$"))


def main():
    res = json.loads((MF / "test_h3_mfbml.json").read_text())["results"]
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.9))
    c10 = fs.level_colour(10)
    for ax, (q, mk, stat, name, unit) in zip(axes, PANELS):
        p = res[q]["points"]
        obs = np.array([x["obs"] for x in p])
        mf = np.array([x["mf_pred"] for x in p])
        err = np.array([mf - [x["mf_lo"] for x in p], [x["mf_hi"] for x in p] - mf])
        kw = {k: v for k, v in fs.series_kw(c10, mk, stat=stat, ours=True).items() if k != "ls"}
        ax.errorbar(obs, mf, yerr=err, ls="none", ecolor=c10, elinewidth=0.8, capsize=2, **kw)
        lo, hi = obs.min() / 1.5, obs.max() * 1.5
        xs = np.geomspace(lo, hi, 20)
        ax.plot(xs, xs, color=fs.NEUTRAL, lw=1.0)
        for f in (0.9, 1.1):
            ax.plot(xs, f * xs, color=fs.NEUTRAL, lw=0.8, ls="--")
        ax.set(xscale="log", yscale="log", xlim=(lo, hi), ylim=(lo, hi), aspect="equal")
        ax.set_xlabel(f"{name} ({unit}), L10 simulation")
        ax.set_ylabel(f"{name} ({unit}), prediction")
        s = res[q]["stats"]
        ax.set_title(f"{name}: error {100 * s['mf_rmse']:.1f}%, {s['mf_cover']}/{s['n']} in 95% band", fontsize=9)
        ax.grid(**fs.GRID_KW)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = MF / f"fig_h3_parity.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print("wrote", out)


if __name__ == "__main__":
    main()
