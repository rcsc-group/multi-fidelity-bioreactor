"""Slide 8 ("result 2"): held-out error against training cost (core-h), log axis.

Data read directly (mtime at writing, `stat -c '%y %n'`):
  /oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/experiments/multifidelity/h3_baselines_learning_curve.json
      (2026-10-21 09:39; learning curve over random subsets of the 13 L10 runs, 15 draws each; tau_95 and EDR)
      NOTE: this file predates the 2026-10-22 switch of EDR to Kim's statistic (edr_max_t). The EDR curve here is the
      time mean (edr_mean_t, which the 2026-10-23 entry still counts as valid). tau_95 is unchanged.
  /oscar/scratch/eaguerov/tmp/cost_audit/cost_table.tsv  (2026-11-04 09:22; sacct)
      cold L8 median (n = 232), cold L9 median (n = 37), median of the 13 warm L10 training runs (names from
      scripts/test_h1_mfbml.hf_points). Cold L10: one run only (4,522 core-h).
Cost model (core-h of the runs each method needs for one fit):
  chain      = 40 x C8 + 16 x C9 + n x C10     (L8 grid, L9 grid, n L10 runs)
  L9 x ratio = 16 x C9 + n x C10              (L9 grid, n L10 runs for the ratio)
  L10 only   = n x C10
N8 = 40 L8 runs and N9 = 16 L9 runs (diary 2026-10-22, periodicity check). Output: figs/s8_cost.png
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from scripts.figstyle import level_colour, MK, series_kw, rcparams, NEUTRAL  # noqa: E402
import scripts.test_h1_mfbml as H1  # noqa: E402

OUT = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026/figs/s8_cost.png")
OUT.parent.mkdir(parents=True, exist_ok=True)
N8, N9 = 40, 16

cost = {}
with open("/oscar/scratch/eaguerov/tmp/cost_audit/cost_table.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        cost.setdefault(r["run"], []).append(r)
cold = {8: [], 9: [], 10: []}
with open("/oscar/scratch/eaguerov/tmp/cost_audit/cost_table.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["cold"] == "True" and int(r["fid"]) in cold:
            cold[int(r["fid"])].append(float(r["corehc_all"]))
C8 = float(np.median(cold[8]))
C9 = float(np.median(cold[9]))
train_names = [H1.hf_run(r, t) for r, t in H1.hf_points()]
warm = []
for nm in train_names:
    rows = cost.get(nm)
    if not rows:
        print("  no cost-table row for", nm)
        continue
    warm.append(sum(float(x["corehc_all"]) for x in rows))
C10 = float(np.median(warm))
print(f"C8 (cold median, n={len(cold[8])}) = {C8:.1f}; C9 (cold median, n={len(cold[9])}) = {C9:.1f}")
print(f"C10 (median of {len(warm)} warm L10 training runs) = {C10:.1f}; values {np.round(sorted(warm), 0)}")

lc = json.load(open(REPO / "experiments/multifidelity/h3_baselines_learning_curve.json"))
print("json keys:", list(lc.keys()))
qk = {"tau_95": "tau_95_t", "edr": "edr_mean_t"}
for k in qk.values():
    if k not in lc:
        raise SystemExit(f"missing key {k}; keys: {list(lc.keys())}")

ns = sorted(int(n) for n in lc["tau_95_t"]["learning_curve"])
print("n:", ns)
cost_fn = {
    "chain": lambda n: N8 * C8 + N9 * C9 + n * C10,
    "l9": lambda n: N9 * C9 + n * C10,
    "l10": lambda n: n * C10,
}
plt.rcParams.update(rcparams())
plt.rcParams.update({"font.size": 14})
COL = {"chain": level_colour(10), "l9": level_colour(9), "l10": NEUTRAL}
LAB = {"chain": "chain (L8 to L9 to L10)", "l9": r"L9 $\times$ ratio", "l10": "L10 only"}
KEY = {"chain": "chain", "l9": "l9_const", "l10": "l10_only"}

fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9), dpi=200, sharey=True)
for ax, (qname, key, title, mk) in zip(axes, [("tau_95", "tau_95_t", r"$\tau_{95}$", MK["tau"]), ("edr", "edr_mean_t", r"$\varepsilon$ time mean", MK["ediss"])]):
    for m in ("chain", "l9", "l10"):
        xs = np.array([cost_fn[m](n) for n in ns])
        ys = np.array([100 * lc[key]["learning_curve"][str(n)][KEY[m]] for n in ns])
        ax.plot(xs, ys, **series_kw(COL[m], mk, stat="mean", ours=True, ms=6.5))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1000, 25000)
    ax.set_xticks([2000, 5000, 10000, 20000])
    ax.set_xticklabels(["2,000", "5,000", "10,000", "20,000"])
    ax.minorticks_off()
    ax.set_title(title, fontsize=16, pad=4)
    ax.set_xlabel("training cost (core-h)", fontsize=14)
    ax.grid(which="major", ls=":", alpha=0.4)
    ax.tick_params(labelsize=13)
axes[0].set_ylabel("error (relative RMS, %)", fontsize=14)
handles = [Line2D([], [], **{**series_kw(COL[m], MK["tau"], stat="mean", ours=True, ms=6.5)}) for m in ("chain", "l9", "l10")]
fig.legend(handles, [LAB["chain"], LAB["l9"], LAB["l10"]], loc="lower center", ncol=3, frameon=False, fontsize=14, bbox_to_anchor=(0.5, -0.03))
fig.subplots_adjust(left=0.1, right=0.98, top=0.9, bottom=0.27, wspace=0.12)
fig.savefig(OUT, dpi=200)
print("wrote", OUT)
for n in ns:
    print(n, {m: round(cost_fn[m](n)) for m in cost_fn}, "tau95 chain/l9/l10 %:",
          [round(100 * lc["tau_95_t"]["learning_curve"][str(n)][KEY[m]], 1) for m in ("chain", "l9", "l10")],
          "edr_mean chain/l9/l10 %:", [round(100 * lc["edr_mean_t"]["learning_curve"][str(n)][KEY[m]], 1) for m in ("chain", "l9", "l10")])
