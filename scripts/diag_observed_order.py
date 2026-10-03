"""Observed order of grid convergence of dtmix, per rpm, for three output transforms.

For three consecutive levels with values z1, z2, z3 (h halves each time), the
increment ratio R = (z2 - z1) / (z3 - z2) is 2^p if z = z0 + a h^p.
  R > 1      monotone convergence, observed order p = log2 R
  0 < R <= 1 increments do not shrink: no convergence on these levels
  R < 0      oscillatory
Transforms: dtmix itself, log dtmix, rate 1/dtmix. Triplets L6-7-8 and L7-8-9.
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs                       # noqa: E402
from scripts.pilot_truncated_chi import RPMS             # noqa: E402

FAM = {6: "kmix_l6_rpm{:g}", 7: "kmix_l7_rpm{:g}", 8: "kmix_l8_rpm{:g}", 9: "fig9_l9_rpm{:g}"}
QOI = sys.argv[1] if len(sys.argv) > 1 else "dtmix"
if QOI == "dtmix":
    TR = {"dtmix": lambda y: y, "log dtmix": np.log, "1/dtmix": lambda y: 1.0 / y}
    CHIS = (0.50, 0.75, 0.95)
    KEY = "dtmix_{:.2f}"
else:                                   # kLa, whole-period estimator, C* = 10/25/50%
    TR = {"kLa": lambda y: y, "log kLa": np.log, "1/kLa": lambda y: 1.0 / y}
    CHIS = (10, 25, 50)
    KEY = "kLa_1T_{}"


def value(lv, r, chi):
    f = ROOT / "runs" / FAM[lv].format(r) / "results.json"
    return json.load(open(f)).get(KEY.format(chi), np.nan) if f.exists() else np.nan


def main():
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(len(TR), len(CHIS), figsize=(12, 7.5), sharex=True, sharey=True)
    for j, chi in enumerate(CHIS):
        for i, (name, g) in enumerate(TR.items()):
            ax = axes[i, j]
            for trip, mk in (((6, 7, 8), "o"), ((7, 8, 9), "s")):
                R = []
                for r in RPMS:
                    z = g(np.array([value(l, r, chi) for l in trip]))
                    R.append((z[1] - z[0]) / (z[2] - z[1]))
                R = np.array(R)
                ok = np.isfinite(R)
                print(f"chi={chi} {name:9s} L{trip[0]}-{trip[2]}: median R = {np.nanmedian(R):6.2f}, "
                      f"monotone-convergent (R>1) {np.mean(R[ok] > 1):.0%}, "
                      f"median p_obs (R>1 only) = {np.median(np.log2(R[ok & (R > 1)])) if (R[ok] > 1).any() else np.nan:.2f}")
                ax.plot(RPMS, R, mk, mfc="w",
                        label=f"L{trip[0]}-L{trip[2]}")
            ax.axhline(1.0, color="0.4", lw=0.8, ls="--")           # R = 1: no convergence
            ax.axhline(2 ** 1.5, color="#d62728", lw=0.8, ls=":")   # p = 1.5
            ax.axhline(0.0, color="0.7", lw=0.6)
            ax.set_yscale("symlog", linthresh=1.0)
            ax.set_yticks([-10, -1, 0, 1, 2.83, 10], ["-10", "-1", "0", "1", "2.8", "10"])
            ax.minorticks_off()
            ax.grid(**fs.GRID_KW)
            if i == 0:
                ax.set_title((rf"$\chi = {chi}$" if QOI == "dtmix" else rf"$C^* = {chi}\%$"), fontsize=10, loc="left")
            if j == 0:
                ax.set_ylabel(f"{name}\nincrement ratio R")
    for ax in axes[-1]:
        ax.set_xlabel("rocking frequency (rpm)")
    axes[0, -1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    out = ROOT / f"experiments/multifidelity/observed_order_{QOI}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
