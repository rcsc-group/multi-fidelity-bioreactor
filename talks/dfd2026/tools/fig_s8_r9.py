"""Slide 8, round 9: x axis in L10-run equivalents; warm-L9 accounting (round 9 addendum).

x unit: ours = core-h / cost of one measured L10 run of our code (666.7 core-h, own mean, tools/cost_rule_r9.py);
        Kim code = core-h / 120 (its L10 run, assumed), i.e. x = n.
Errors: data/chain_variants_lc.json, data/lc_tau_95_t_edr_max_t.json, data/lc_h1_tau_95_t_edr_max_t.json (unchanged).
Costs: data/cost_rule_r9.json (warm-L9 estimate for the 10 cold fig9 theta-7 runs; data/warm_l9_r9.json).
Outputs: figs/s8_cost.png, data/crossings_r9.json. Prints the table and crossings (core-h and L10-equivalents).
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
CR = json.load(open(TALK / "data" / "cost_rule_r9.json"))
UNIT = CR["L10_unit_core_h"]
NS = (3, 4, 6, 8, 10, 12)
QS = ("tau_95_t", "edr_max_t")
PINK, BLUE, GREY, BLACK = "#CC79A7", "#0072B2", "0.45", "black"
VARIANTS = {"full": "chain_full", "warm9": "chain_warm9", "warm9_L8_20": "chain_warm9_L8_20", "warm9_L8_10": "chain_warm9_L8_10"}


def xval(m, n):
    """Training cost in L10-run equivalents."""
    if m == "Kim_code":
        return n * CR[m]["inc"] / 120.0
    return (CR[m]["base"] + n * CR[m]["inc"]) / UNIT


def core_h(m, n):
    return CR[m]["base"] + n * CR[m]["inc"]


chain_pts = {q: [] for q in QS}
for name, key in VARIANTS.items():
    for q in QS:
        for n in NS:
            chain_pts[q].append(dict(variant=name, n=n, x=xval(key, n), ch=core_h(key, n), err=CV[name]["err"][q][str(n)]))


def envelope(pts):
    out, best = [], np.inf
    for p in sorted(pts, key=lambda p: (p["x"], p["err"])):
        if p["err"] < best:
            out.append(p)
            best = p["err"]
    return out


def first_under(points, thr):
    hits = [p for p in sorted(points, key=lambda p: p["x"]) if p["err"] <= thr]
    return hits[0] if hits else None


plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Inter", "Helvetica Neue", "DejaVu Sans"],
                     "mathtext.fontset": "dejavusans", "font.size": 14, "axes.linewidth": 1.1,
                     "xtick.direction": "in", "ytick.direction": "in"})
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.9), dpi=200, sharey=True)
res = {}
for ax, q, qn in zip(axes, QS, [r"$\tau_{95}$", r"$\varepsilon_{\max}$ (EDR)"]):
    for p in chain_pts[q]:
        ax.plot(p["x"], 100 * p["err"], marker="o", ms=3.5, mfc="none", mec=PINK, alpha=0.35, ls="none")
    env = envelope(chain_pts[q])
    ax.plot([p["x"] for p in env], [100 * p["err"] for p in env], color=PINK, ls="-", lw=1.8, marker="o", ms=5.5, mfc=PINK, mec=PINK)
    full = [(xval("chain_full", n), core_h("chain_full", n), CV["full"]["err"][q][str(n)]) for n in NS]
    ax.plot([x for x, _, _ in full], [100 * e for _, _, e in full], color=PINK, ls=":", lw=1.1, marker="o", ms=5, mfc="white", mec=PINK)
    l10 = [(xval("L10_GP_ours", n), core_h("L10_GP_ours", n), LC[q]["learning_curve"][str(n)]["l10_only"]) for n in NS]
    ax.plot([x for x, _, _ in l10], [100 * e for _, _, e in l10], color=GREY, ls=":", lw=1.4, marker="o", ms=6, mfc="white", mec=GREY)
    kim = [(xval("Kim_code", n), 120.0 * n, LC[q]["learning_curve"][str(n)]["l10_only"]) for n in NS]
    ax.plot([x for x, _, _ in kim], [100 * e for _, _, e in kim], color=BLACK, ls="-", lw=1.4, marker="s", ms=6, mfc=BLACK, mec=BLACK)
    h1 = [(xval("H1", n), core_h("H1", n), LH[q]["learning_curve"][str(n)]["h1"]) for n in NS]
    ax.plot([x for x, _, _ in h1], [100 * e for _, _, e in h1], color=BLUE, ls=":", lw=1.4, marker="o", ms=6, mfc="white", mec=BLUE)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.5, 60)
    ax.set_xticks([1, 2, 5, 10, 20, 50])
    ax.set_xticklabels(["1", "2", "5", "10", "20", "50"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(1.5, 400)
    ax.yaxis.set_major_locator(FixedLocator([2, 3, 5, 10, 20, 50, 100, 200, 400]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_title(qn, fontsize=14, pad=3)
    ax.set_xlabel("training cost (L10-run equivalents)", fontsize=14)
    ax.grid(which="major", ls=":", alpha=0.4)
    ax.tick_params(labelsize=13)
    res[q] = {}
    for thr in (0.10, 0.20):
        cands = {
            "chain envelope": first_under(chain_pts[q], thr),
            "chain, full set": first_under([dict(x=x, err=e, ch=c) for x, c, e in full], thr),
            "L10 GP, ours": first_under([dict(x=x, err=e, ch=c) for x, c, e in l10], thr),
            "L10 GP, Kim code (120/run)": first_under([dict(x=x, err=e, ch=c) for x, c, e in kim], thr),
            "H1 (L8 to L10)": first_under([dict(x=x, err=e, ch=c) for x, c, e in h1], thr),
        }
        res[q][f"{int(thr * 100)}%"] = {k: (None if v is None else dict(x=round(v["x"], 2), core_h=round(v["ch"]) if "ch" in v else None, n=v.get("n"), variant=v.get("variant"))) for k, v in cands.items()}

axes[0].set_ylabel("held-out error (relative RMS, %)", fontsize=14)
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
t10, t20, e10 = res["tau_95_t"]["10%"], res["tau_95_t"]["20%"], res["edr_max_t"]["10%"]
ch = res["tau_95_t"]["10%"]["chain envelope"]
h1_10 = res["tau_95_t"]["10%"]["H1 (L8 to L10)"]
axes[0].text(1.1, 150, f"10 % error, tau 95:\nchain {ch['x']:.1f} L10-eq\nH1 {h1_10['x']:.1f} L10-eq;\nL10 GP not reached" if h1_10 else "", fontsize=11, ha="left", va="center", color="0.1")
ce = res["edr_max_t"]["10%"]["chain, full set"]
axes[1].text(0.8, 9, f"10 % error, EDR:\nonly the chain reaches it\n({ce['x']:.1f} L10-eq)", fontsize=11, ha="left", va="center", color="0.1")
fig.savefig(OUT, dpi=200)
(TALK / "data" / "crossings_r9.json").write_text(json.dumps(res, indent=1))
print("wrote", OUT, "unit", round(UNIT, 1))
for q in res:
    for thr in res[q]:
        print(q, thr)
        for k, v in res[q][thr].items():
            print("   ", f"{k:28s}", v)
