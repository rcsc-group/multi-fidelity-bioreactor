"""Slide 8, round 7: one accounting rule (tools/cost_rule_r7.py) for every curve; envelope and crossings recomputed.

Errors: data/chain_variants_lc.json (the four chain variants, same held-out draws), data/lc_tau_95_t_edr_max_t.json
(approved L10-only GP errors, used for "L10 GP, ours" and for the Kim curve), data/lc_h1_tau_95_t_edr_max_t.json (H1).
Costs: data/cost_rule_r7.json. Cost of n L10 runs = base + n x per-L10 charge (expected value over random subsets).
Kim code: 120 core-h per run, no sources (user's assumption; same errors as L10 GP, ours).
Outputs: figs/s8_cost.png, data/crossings_r7.json. Prints the table and the crossings.
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
OUT = TALK / "figs" / "s8_cost.png"
CV = json.load(open(TALK / "data" / "chain_variants_lc.json"))
LC = json.load(open(TALK / "data" / "lc_tau_95_t_edr_max_t.json"))
LH = json.load(open(TALK / "data" / "lc_h1_tau_95_t_edr_max_t.json"))
CR = json.load(open(TALK / "data" / "cost_rule_r7.json"))
NS = (3, 4, 6, 8, 10, 12)
QS = ("tau_95_t", "edr_max_t")
PINK, BLUE, GREY, BLACK = "#CC79A7", "#0072B2", "0.45", "black"
VARIANTS = {"full": "chain_full", "warm9": "chain_warm9", "warm9_L8_20": "chain_warm9_L8_20", "warm9_L8_10": "chain_warm9_L8_10"}

def cost(m, n):
    return CR[m]["base"] + n * CR[m]["inc"]

chain_pts = {q: [] for q in QS}
for name, key in VARIANTS.items():
    for q in QS:
        for n in NS:
            chain_pts[q].append(dict(variant=name, n=n, cost=cost(key, n), err=CV[name]["err"][q][str(n)]))


def envelope(pts):
    out, best = [], np.inf
    for p in sorted(pts, key=lambda p: (p["cost"], p["err"])):
        if p["err"] < best:
            out.append(p)
            best = p["err"]
    return out


def first_under(points, thr):
    hits = [p for p in sorted(points, key=lambda p: p["cost"]) if p["err"] <= thr]
    return hits[0] if hits else None


res = {}
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14, "axes.linewidth": 1.1,
                     "xtick.direction": "in", "ytick.direction": "in"})
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.9), dpi=200, sharey=True)
for ax, q, qn in zip(axes, QS, [r"$\tau_{95}$", r"$\varepsilon_{\max}$ (EDR)"]):
    for p in chain_pts[q]:
        ax.plot(p["cost"], 100 * p["err"], marker="o", ms=3.5, mfc="none", mec=PINK, alpha=0.35, ls="none")
    env = envelope(chain_pts[q])
    ax.plot([p["cost"] for p in env], [100 * p["err"] for p in env], color=PINK, ls="-", lw=1.8, marker="o", ms=5.5, mfc=PINK, mec=PINK)
    full = [(cost("chain_full", n), CV["full"]["err"][q][str(n)]) for n in NS]
    ax.plot([c for c, _ in full], [100 * e for _, e in full], color=PINK, ls=":", lw=1.1, marker="o", ms=5, mfc="white", mec=PINK)
    l10 = [(cost("L10_GP_ours", n), LC[q]["learning_curve"][str(n)]["l10_only"]) for n in NS]
    ax.plot([c for c, _ in l10], [100 * e for _, e in l10], color=GREY, ls=":", lw=1.4, marker="o", ms=6, mfc="white", mec=GREY)
    kim = [(cost("Kim_code", n), LC[q]["learning_curve"][str(n)]["l10_only"]) for n in NS]
    ax.plot([c for c, _ in kim], [100 * e for _, e in kim], color=BLACK, ls="-", lw=1.4, marker="s", ms=6, mfc=BLACK, mec=BLACK)
    h1 = [(cost("H1", n), LH[q]["learning_curve"][str(n)]["h1"]) for n in NS]
    ax.plot([c for c, _ in h1], [100 * e for _, e in h1], color=BLUE, ls=":", lw=1.4, marker="o", ms=6, mfc="white", mec=BLUE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(300, 30000)
    ax.set_xticks([1000, 5000, 20000])
    ax.set_xticklabels(["1k", "5k", "20k"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(1.5, 400)
    ax.yaxis.set_major_locator(FixedLocator([2, 3, 5, 10, 20, 50, 100, 200, 400]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_title(qn, fontsize=14, pad=3)
    ax.set_xlabel("training cost (core-h)", fontsize=14)
    ax.grid(which="major", ls=":", alpha=0.4)
    ax.tick_params(labelsize=13)
    res[q] = {}
    for thr in (0.10, 0.20):
        cands = {
            "chain envelope": first_under(chain_pts[q], thr),
            "chain, approved set": first_under([dict(cost=c, err=e) for c, e in full], thr),
            "L10 GP, ours": first_under([dict(cost=c, err=e) for c, e in l10], thr),
            "L10 GP, Kim code (120/run)": first_under([dict(cost=c, err=e) for c, e in kim], thr),
            "H1 (L8 to L10)": first_under([dict(cost=c, err=e) for c, e in h1], thr),
        }
        res[q][f"{int(thr * 100)}%"] = {k: (None if v is None else dict(cost=round(v["cost"]), n=v.get("n"), variant=v.get("variant"))) for k, v in cands.items()}

axes[0].set_ylabel("error (relative RMS, %)", fontsize=14)
axes[0].text(1250, 150, "10 % error, tau 95:\nchain 4,582 core-h\nH1 16,700 core-h;\nL10 GP not reached", fontsize=11, ha="left", va="center", color="0.1")
axes[1].text(700, 9, "10 % error, EDR:\nonly the chain reaches it\n(13,246 core-h)", fontsize=11, ha="left", va="center", color="0.1")
handles = [
    Line2D([], [], color=PINK, ls="-", lw=1.8, marker="o", mfc=PINK, mec=PINK, ms=5.5),
    Line2D([], [], color=PINK, ls="none", marker="o", mfc="none", mec=PINK, alpha=0.5, ms=4),
    Line2D([], [], color=GREY, ls=":", lw=1.4, marker="o", mfc="white", mec=GREY, ms=6),
    Line2D([], [], color=BLACK, ls="-", lw=1.4, marker="s", mfc=BLACK, mec=BLACK, ms=6),
    Line2D([], [], color=BLUE, ls=":", lw=1.4, marker="o", mfc="white", mec=BLUE, ms=6),
]
labels = ["chain, lower envelope of variants", "chain variants (faint)", "L10 GP, ours", "L10 GP, Kim code (assumed 120 core-h/run)", "H1, L8 to L10"]
fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, fontsize=12, bbox_to_anchor=(0.5, 1.0), handlelength=2.0, columnspacing=1.0)
fig.subplots_adjust(left=0.11, right=0.985, top=0.80, bottom=0.14, wspace=0.12)
fig.savefig(OUT, dpi=200)
(TALK / "data" / "crossings_r7.json").write_text(json.dumps(res, indent=1))
print("wrote", OUT)
print("per-method cost (base, per L10 run):")
for k in ("chain_full", "chain_warm9", "chain_warm9_L8_20", "chain_warm9_L8_10", "L10_GP_ours", "H1", "Kim_code"):
    print(f"   {k:20s} base {CR[k]['base']:9.1f}  per L10 {CR[k]['inc']:8.1f}")
for q in res:
    for thr in res[q]:
        print(q, thr)
        for k, v in res[q][thr].items():
            print("   ", f"{k:28s}", v)
