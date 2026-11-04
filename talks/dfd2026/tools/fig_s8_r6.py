"""Slide 8, round 6 (A2-A4): chain lower envelope over cheaper variants, Kim-code L10 GP (assumed cost), crossings.

Inputs (mtime from `stat -c '%y %n'`):
  data/chain_variants_lc.json      (2026-11-04 11:2x; tools/chain_variants_lc.py: 4 chain variants, same subsets)
  data/lc_tau_95_t_edr_max_t.json  (2026-11-04 10:20; approved L10-only GP curve, key l10_only; chain, approved)
  data/lc_h1_tau_95_t_edr_max_t.json (2026-11-04 10:4x; H1)
  data/costs.json                  (sacct per run id; L10-only convention 1,244 core-h per run incl. source)
Kim: L10 GP error values are the same as ours; cost = n x 120 core-h (ASSUMED, the user's value; KIM_CORE_H_PER_CASE 120 is
  marked PENDING AUDIT in scripts/campaign2d/campaign.py).
Output: figs/s8_cost.png, data/crossings_r6.json
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
CS = json.load(open(TALK / "data" / "costs.json"))
NS = (3, 4, 6, 8, 10, 12)
QS = ("tau_95_t", "edr_max_t")
KIM_PER_RUN = 120.0
L10_GP_PER_RUN = CS["inc_l10"]
PINK, BLUE, GREY, BLACK = "#CC79A7", "#0072B2", "0.45", "black"

# chain points: (cost, error) for every variant and n
chain_pts = {q: [] for q in QS}
for name, v in CV.items():
    for q in QS:
        for n in NS:
            cost = v["costs"]["base"] + n * v["costs"]["inc_per_L10"]
            chain_pts[q].append(dict(variant=name, n=n, cost=cost, err=v["err"][q][str(n)]))


def frontier(pts):
    """Lower envelope: points not dominated by a cheaper point with lower error."""
    out, best = [], np.inf
    for p in sorted(pts, key=lambda p: (p["cost"], p["err"])):
        if p["err"] < best:
            out.append(p)
            best = p["err"]
    return out


def curve(q, key_src=None):
    return [(n, LC[q]["learning_curve"][str(n)]) for n in NS]


def crossings(points, thr):
    """First point (sorted by cost) with error <= thr. points: list of (cost, err, label)."""
    hits = [p for p in sorted(points, key=lambda p: p[0]) if p[1] <= thr]
    return (round(hits[0][0]), hits[0][2]) if hits else None


res = {}
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.9), dpi=200, sharey=True)
for ax, q, qn in zip(axes, QS, [r"$\tau_{95}$", r"$\varepsilon_{\max}$ (EDR)"]):
    # chain variants (faint) and envelope (pink)
    for p in chain_pts[q]:
        ax.plot(p["cost"], 100 * p["err"], marker="o", ms=3.5, mfc="none", mec=PINK, alpha=0.35, ls="none")
    env = frontier(chain_pts[q])
    ax.plot([p["cost"] for p in env], [100 * p["err"] for p in env], color=PINK, ls="-", lw=1.8, marker="o", ms=5.5, mfc=PINK, mec=PINK)
    # approved chain (full set), for reference
    full = [(CV["full"]["costs"]["base"] + n * CV["full"]["costs"]["inc_per_L10"], CV["full"]["err"][q][str(n)]) for n in NS]
    ax.plot([c for c, _ in full], [100 * e for _, e in full], color=PINK, ls=":", lw=1.1, marker="o", ms=5, mfc="white", mec=PINK)
    # ours: L10-only GP (approved values, cost 1,244 per run incl. source)
    l10 = [(n * L10_GP_PER_RUN, LC[q]["learning_curve"][str(n)]["l10_only"]) for n in NS]
    ax.plot([c for c, _ in l10], [100 * e for _, e in l10], color=GREY, ls=":", lw=1.4, marker="o", ms=6, mfc="white", mec=GREY)
    # Kim code: same error values, cost 120 per run (assumed)
    kim = [(n * KIM_PER_RUN, LC[q]["learning_curve"][str(n)]["l10_only"]) for n in NS]
    ax.plot([c for c, _ in kim], [100 * e for _, e in kim], color=BLACK, ls="-", lw=1.4, marker="s", ms=6, mfc=BLACK, mec=BLACK)
    # H1
    h1 = [(CS["base_h1"] + n * CS["inc_h1"], LH[q]["learning_curve"][str(n)]["h1"]) for n in NS]
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

    # crossings at 10 % and 20 % (grid points only; "not reached" if no n in NS qualifies)
    for thr in (0.10, 0.20):
        ent = {
            "chain envelope": crossings([(p["cost"], p["err"], f"{p['variant']} n={p['n']}") for p in chain_pts[q]], thr),
            "chain approved (40 L8, 16 L9)": crossings([(c, e, "full") for c, e in full], thr),
            "L10 GP, ours (1,244/run)": crossings([(c, e, "l10") for c, e in l10], thr),
            "L10 GP, Kim code (120/run, assumed)": crossings([(c, e, "kim") for c, e in kim], thr),
            "H1 (L8 to L10)": crossings([(c, e, "h1") for c, e in h1], thr),
        }
        res.setdefault(q, {})[f"{int(thr * 100)}%"] = ent

axes[0].set_ylabel("error (relative RMS, %)", fontsize=14)
axes[0].text(1700, 110, "same 20 % error (tau 95):\nchain envelope 2,918 core-h\nKim code, assumed 1,200 core-h", fontsize=11, ha="left", va="center", color="0.1")
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
(TALK / "data" / "crossings_r6.json").write_text(json.dumps(res, indent=1))
print("wrote", OUT)
for q in res:
    for thr in res[q]:
        print(q, thr)
        for k, v in res[q][thr].items():
            print("   ", k, v)
