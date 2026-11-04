"""Slide 8 ("result 2"), round 2: held-out error against ACTUAL training cost (sacct, per run id), log axis.

Inputs (mtime from `stat -c '%y %n'`):
  data/lc_tau_95_t_edr_max_t.json  (2026-11-04 10:20; made by tools/h3_lc_copy.py with the repo's learning-curve code,
      seed 20261021, 15 draws per n; the tau_95 values reproduce experiments/multifidelity/h3_baselines_learning_curve.json
      (2026-10-21 09:39) exactly; the EDR column is Kim's edr_max_t, the statistic used on slide 7)
  data/costs.json  (2026-11-04 10:2x; tools/cost_s8.py: sum of corehc_all per run id, cost_table.tsv 2026-11-04 09:22)

Cost model (core-h of the runs each training set needs; each run id counted once, sources of warm starts included):
  chain       = all 40 L8 + all 16 L9 (+ their sources) + n L10 runs (own cost; their sources are already in the base)
  L9 x ratio  = all 16 L9 (+ their L8 sources) + n L10 runs (own cost)
  L9 x log-linear ratio = same cost as L9 x ratio (it needs the same runs; n >= 4 only)
  L10 only    = n L10 runs, each with its own warm-start source (L9 fig9 run or L8 t6 run)

Outputs: figs/s8_cost.png and data/crossings.json (cost at which each method first reaches 10 % and 5 % error).
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
OUT = TALK / "figs" / "s8_cost.png"
lc = json.load(open(TALK / "data" / "lc_tau_95_t_edr_max_t.json"))
cs = json.load(open(TALK / "data" / "costs.json"))
base_chain, inc_chain = cs["base_chain"], cs["inc_chain"]
base_ratio, inc_ratio = cs["base_ratio"], cs["inc_ratio"]
inc_l10 = cs["inc_l10"]
COST = {
    "chain": lambda n: base_chain + n * inc_chain,
    "ratio": lambda n: base_ratio + n * inc_ratio,
    "loglin": lambda n: base_ratio + n * inc_ratio,
    "l10": lambda n: n * inc_l10,
}
KEY = {"chain": "chain", "ratio": "l9_const", "loglin": "l9_loglin", "l10": "l10_only"}
COL = {"chain": "#CC79A7", "ratio": "#009E73", "loglin": "#009E73", "l10": "0.45"}
LS = {"chain": ":", "ratio": ":", "loglin": "-.", "l10": ":"}
MK = {"tau_95_t": "o", "edr_max_t": "s"}
NAME = {"chain": "chain (L8, L9, L10)", "ratio": "L9 x constant ratio", "loglin": "L9 x log-linear ratio", "l10": "L10 only"}
QN = {"tau_95_t": r"$\tau_{95}$", "edr_max_t": r"$\varepsilon_{\max}$ (EDR)"}

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14, "axes.linewidth": 1.1,
                     "xtick.direction": "in", "ytick.direction": "in"})

ns = sorted(int(n) for n in lc["tau_95_t"]["learning_curve"])
cross = {}
fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6), dpi=200, sharey=True)
for ax, q in zip(axes, ["tau_95_t", "edr_max_t"]):
    cross[q] = {}
    for m in ("chain", "ratio", "loglin", "l10"):
        pts = [(n, lc[q]["learning_curve"][str(n)].get(KEY[m])) for n in ns]
        pts = [(n, e) for n, e in pts if e is not None]
        xs = np.array([COST[m](n) for n, _ in pts])
        ys = np.array([100 * e for _, e in pts])
        ax.plot(xs, ys, color=COL[m], ls=LS[m], lw=1.4, marker=MK[q], ms=6.5, mfc="white", mec=COL[m], mew=1.3)
        cross[q][m] = {}
        for thr in (10, 5):
            hit = [(n, COST[m](n)) for n, e in pts if 100 * e <= thr]
            cross[q][m][f"{thr}%"] = (dict(n=hit[0][0], cost_core_h=round(hit[0][1])) if hit else "not reached by n = 12")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1000, 30000)
    ax.set_xticks([2000, 5000, 10000, 20000])
    ax.set_xticklabels(["2,000", "5,000", "10,000", "20,000"])
    ax.minorticks_off()
    ax.set_title(QN[q], fontsize=14, pad=3)
    ax.set_xlabel("training cost (core-h)", fontsize=14)
    ax.grid(which="major", ls=":", alpha=0.4)
    ax.tick_params(labelsize=13)
axes[0].set_ylabel("error (relative RMS, %)", fontsize=14)
axes[0].set_yticks([2, 5, 10, 20, 50, 100, 200])
axes[0].set_yticklabels(["2", "5", "10", "20", "50", "100", "200"])
axes[1].set_yticks([2, 5, 10, 20, 50, 100, 200])
axes[1].set_yticklabels([])
handles = [Line2D([], [], color=COL[m], ls=LS[m], lw=1.4, marker="o", mfc="white", mec=COL[m], ms=6.5) for m in ("chain", "ratio", "loglin", "l10")]
fig.legend(handles, [NAME[m] for m in ("chain", "ratio", "loglin", "l10")], loc="upper center", ncol=4, frameon=False,
           fontsize=13, bbox_to_anchor=(0.5, 1.01), handlelength=2.2, columnspacing=1.0)
fig.subplots_adjust(left=0.1, right=0.985, top=0.8, bottom=0.2, wspace=0.12)
fig.savefig(OUT, dpi=200)
json.dump(dict(cross=cross, base_chain=base_chain, inc_chain=inc_chain, base_ratio=base_ratio, inc_ratio=inc_ratio, inc_l10=inc_l10),
          open(TALK / "data" / "crossings.json", "w"), indent=1)
print("wrote", OUT)
for q in cross:
    for m in cross[q]:
        print(q, m, cross[q][m])
