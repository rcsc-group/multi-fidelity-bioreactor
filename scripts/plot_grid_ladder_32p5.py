"""Grid ladder at 32.5 rpm, 7 deg: Kim's quantities vs level, with the E&H limit.

Each panel: our L6-L9 values (markers), the Eca & Hoekstra (2014) extrapolated
value phi0 at h = 0 with L9's 95% bound U drawn on the L9 point, and Kim's value
(solid line). x is the cell size relative to L10, h = 2^(10-N).
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs                       # noqa: E402
from scripts.eca_hoekstra import numerical_uncertainty  # noqa: E402

RUNS = {9: "fig9_l9_rpm32.5", 8: "fig9_ladder_l8_rpm32.5",
        7: "fig9_ladder_l7_rpm32.5", 6: "fig9_validate_l6_rpm32.5"}
KIM = pd.read_csv(ROOT / "experiments/kimetal2024/csv_raw/mixing_kla_vs_frequency.csv").set_index("RPM").loc[32.5]
PANELS = [("kLa_1T_25", "kLa_exp5pts_25", r"$k_La$ at $C^*=25\%$ (h$^{-1}$)"),
          ("kLa_1T_50", "kLa_exp5pts_50", r"$k_La$ at $C^*=50\%$ (h$^{-1}$)"),
          ("dtmix_0.50", "dtmix_strict_0.5", r"$\Delta t_{\mathrm{mix}}$ at $\chi=0.50$ (s)"),
          ("dtmix_0.95", "dtmix_strict_0.95", r"$\Delta t_{\mathrm{mix}}$ at $\chi=0.95$ (s)")]

plt.rcParams.update(fs.rcparams())
fig, axes = plt.subplots(1, 4, figsize=(13.0, 3.3))
for ax, (key, kcol, lab) in zip(axes, PANELS):
    lv = [9, 8, 7, 6]
    q = np.array([json.load(open(ROOT / "runs" / RUNS[N] / "results.json"))[key] for N in lv])
    h = np.array([2.0 ** (10 - N) for N in lv])
    r = numerical_uncertainty(h, q, index=0)
    for N, hh, v in zip(lv, h, q):
        ax.plot(hh, v, "o", color=fs.level_colour(N), ms=6, zorder=3)
    ax.errorbar(h[0], q[0], yerr=[[min(r["U"], 0.999 * q[0])], [r["U"]]], fmt="none",
                ecolor=fs.level_colour(9), capsize=3, lw=1.1)
    hg = np.linspace(0, h.max(), 100)
    ax.plot(hg, r["fit"](hg), color="0.5", lw=1.0, ls="--")
    if r["phi0"] > 0:
        ax.plot(0, r["phi0"], "D", mfc="w", mec="0.3", ms=6, zorder=3)
    ax.axhline(float(KIM[kcol]), color=fs.KIM, lw=fs.LW_KIM)
    ax.set_xlabel(r"cell size $h/h_{10}$")
    ax.set_title(lab, loc="left", fontsize=9.5)
    ax.set_xscale("symlog", linthresh=1.0)
    ax.set_xticks([0, 1, 2, 4, 8, 16])
    ax.set_xticklabels(["0", "1", "2", "4", "8", "16"])
    ax.grid(**fs.GRID_KW)
from matplotlib.lines import Line2D
handles = [Line2D([], [], ls="none", marker="o", color=fs.level_colour(N), label=f"L{N}") for N in (6, 7, 8, 9)]
handles += [Line2D([], [], color="0.5", ls="--", label="E&H fit"),
            Line2D([], [], ls="none", marker="D", mfc="w", mec="0.3", label=r"extrapolated, $h\to0$"),
            Line2D([], [], color=fs.KIM, lw=fs.LW_KIM, label="Kim et al.")]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 0.92), frameon=False, fontsize=8)
fig.suptitle(r"$f_b=32.5$ rpm, $\theta_{b,\max}=7^\circ$", x=0.01, ha="left", fontsize=10)
fig.tight_layout()
out = ROOT / "experiments/multifidelity/grid_ladder_32p5.png"
fig.savefig(out, dpi=170, bbox_inches="tight")
print("saved", out)
