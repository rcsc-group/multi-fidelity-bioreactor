"""Slide 8 ("result 2"), round 3: held-out error against actual training cost, three methods, log cost axis.

Methods: chain (L8 -> L9 -> L10), L10 only (single-fidelity GP), H1 (two-level, direct L8 -> L10, raw scale).
The L9 ratio curves are removed (lead decision, round 3).

Inputs (mtime from `stat -c '%y %n'`):
  data/lc_tau_95_t_edr_max_t.json  (2026-11-04 10:20; chain and L10-only; tools/h3_lc_copy.py; seed 20261021; DRAWS 15; n = 3..12)
  data/lc_h1_tau_95_t_edr_max_t.json (2026-11-04 10:4x; tools/h1_lc_copy.py; same seed, same draws, raw scale)
  data/costs.json (2026-11-04 10:4x; tools/cost_s8.py; sum of corehc_all per run id from cost_table.tsv 2026-11-04 09:22)

Cost model (core-h, each run id counted once, warm-start sources included):
  chain : 40 L8 + 16 L9 (+ sources) base, then + n x 667 (own L10 cost; sources already in base)
  L10 only : n x 1,244 (each L10 run with its own source: L9 fig9 run or L8 t6 run)
  H1    : 40 L8 base (2,043), then + n x 1,222 (each L10 run with its source, unless the source is an L8 run already in base)
Outputs: figs/s8_cost.png, data/crossings_r3.json (cost at which each method first reaches 10 % and 5 % error)
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, FixedLocator, NullFormatter

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
OUT = TALK / "figs" / "s8_cost.png"
lc = json.load(open(TALK / "data" / "lc_tau_95_t_edr_max_t.json"))
lh1 = json.load(open(TALK / "data" / "lc_h1_tau_95_t_edr_max_t.json"))
cs = json.load(open(TALK / "data" / "costs.json"))
COST = {
    "chain": lambda n: cs["base_chain"] + n * cs["inc_chain"],
    "l10": lambda n: n * cs["inc_l10"],
    "h1": lambda n: cs["base_h1"] + n * cs["inc_h1"],
}
KEY = {"chain": ("chain", lc), "l10": ("l10_only", lc), "h1": ("h1", lh1)}
COL = {"chain": "#CC79A7", "l10": "0.45", "h1": "#0072B2"}
MK = {"tau_95_t": "o", "edr_max_t": "s"}
NAME = {"chain": "chain (L8, L9, L10)", "l10": "L10 only (GP)", "h1": "two-level, L8 to L10 (H1)"}
QN = {"tau_95_t": r"$\tau_{95}$", "edr_max_t": r"$\varepsilon_{\max}$ (EDR)"}

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14, "axes.linewidth": 1.1,
                     "xtick.direction": "in", "ytick.direction": "in"})
ns = sorted(int(n) for n in lc["tau_95_t"]["learning_curve"])
cross = {}
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.9), dpi=200, sharey=True)
for ax, q in zip(axes, ["tau_95_t", "edr_max_t"]):
    cross[q] = {}
    for m in ("chain", "l10", "h1"):
        key, src = KEY[m]
        pts = [(n, src[q]["learning_curve"][str(n)][key]) for n in ns if key in src[q]["learning_curve"][str(n)]]
        xs = np.array([COST[m](n) for n, _ in pts])
        ys = np.array([100 * e for _, e in pts])
        ax.plot(xs, ys, color=COL[m], ls=":", lw=1.5, marker=MK[q], ms=6.5, mfc="white", mec=COL[m], mew=1.3)
        cross[q][m] = {}
        for thr in (10, 5):
            hit = [(n, COST[m](n)) for n, e in pts if 100 * e <= thr]
            cross[q][m][f"{thr}%"] = (dict(n=hit[0][0], cost_core_h=round(hit[0][1])) if hit else "not reached by n = 12")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1000, 40000)
    ax.set_xticks([2000, 5000, 10000, 20000])
    ax.set_xticklabels(["2,000", "5,000", "10,000", "20,000"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(3, 400)
    ax.yaxis.set_major_locator(FixedLocator([5, 10, 20, 50, 100, 200, 400]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_title(QN[q], fontsize=14, pad=3)
    ax.set_xlabel("training cost (core-h)", fontsize=14)
    ax.grid(which="major", ls=":", alpha=0.4)
    ax.tick_params(labelsize=13)
axes[0].set_ylabel("error (relative RMS, %)", fontsize=14)
handles = [Line2D([], [], color=COL[m], ls=":", lw=1.5, marker="o", mfc="white", mec=COL[m], ms=6.5) for m in ("chain", "l10", "h1")]
fig.legend(handles, [NAME[m] for m in ("chain", "l10", "h1")], loc="upper center", ncol=3, frameon=False,
           fontsize=13, bbox_to_anchor=(0.5, 1.0), handlelength=2.2, columnspacing=1.2)
fig.subplots_adjust(left=0.11, right=0.985, top=0.84, bottom=0.14, wspace=0.12)
fig.savefig(OUT, dpi=200)
json.dump(dict(cross=cross, base_chain=cs["base_chain"], inc_chain=cs["inc_chain"], inc_l10=cs["inc_l10"],
               base_h1=cs["base_h1"], inc_h1=cs["inc_h1"]), open(TALK / "data" / "crossings_r3.json", "w"), indent=1)
print("wrote", OUT)
for q in cross:
    for m in cross[q]:
        print(q, m, cross[q][m])
