"""Two figures for the 2026-11-10 discussion.
(1) levels_theta7.png: kLa_1T_10 and t_mix,25 vs rpm at theta = 7 (the 1D line), L7 / L8 / L9, every real run in the
    frozen snapshots. Grammar: colour = level; marker = quantity (circle kLa, triangle t_mix); dotted = ours.
(2) err_vs_coreh.png: held-out L9 RMSE (median of 20 site draws, mf_ablation.csv) vs training core-h, per quantity.
    Cost charges EVERY run the model uses (measured sacct core-h): L9-only n9 * (200 + 16.8 spin-up);
    chain n9 * 200 + (all L8 rows) * 18 + (all L7 rows, 3-level only) * 5. Colour = lowest level used.
  uv run python scripts/bo2d/fig_levels_theta7_and_cost.py OUTDIR"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
S = ROOT / "scripts/bo2d/state"
plt.rcParams.update(fs.rcparams())
snap = {q: pd.read_csv(S / f"truth_snapshot_{q}.csv") for q in ("kLa_1T_10", "dtmix_0.25", "tau_95")}
for d in snap.values():
    d.drop(d[~(np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0))].index, inplace=True)

fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.2))
for ax, (q, mk, lab) in zip(axs, (("kLa_1T_10", "o", "kLa (1/h)"), ("dtmix_0.25", "^", "t_mix,25 (s)"))):
    d = snap[q]; d = d[np.isclose(d.theta, 7.0, atol=0.05)]
    for lev in (7, 8, 9):
        s = d[d.level == lev].groupby("rpm").y.mean().sort_index()
        ax.plot(s.index, s.values, ls=":", marker=mk, color=fs.LEVEL_COLOUR[lev], lw=1.3, ms=5, label=f"L{lev}")
    ax.set_yscale("log"); ax.set_xlabel("rpm"); ax.set_ylabel(lab)
axs[1].legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
fig.tight_layout(); fig.savefig(out / "levels_theta7.png", dpi=200, bbox_inches="tight")

a = pd.read_csv(S / "mf_ablation.csv")
n8 = {q: int((snap[q].level == 8).sum()) for q in snap}; n7 = {q: int((snap[q].level == 7).sum()) for q in snap}


def cost(v, q, n9):
    if v in ("U1", "B1"):
        return n9 * 216.8
    lv = 3 if v[1] == "3" else 2
    return n9 * 200 + n8[q] * 18 + (n7[q] * 5 if lv == 3 else 0)


STY = {"B1": (9, "o"), "B2-lin": (8, "o"), "B3-lin": (7, "o")}
TITLE = {"kLa_1T_10": "kLa", "dtmix_0.25": "mixing time", "tau_95": "shear stress"}
fig, axs = plt.subplots(1, 3, figsize=(10, 3.2))
for ax, q in zip(axs, TITLE):
    g = a[(a.qoi == q) & a.variant.isin(STY)].groupby(["variant", "n9"]).rmse.median()
    for v, (lev, mk) in STY.items():
        s = g[v]; x = [cost(v, q, n) for n in s.index]
        ax.plot(x, s.values, ls="-" if v == "B1" else ":", marker=mk, color=fs.LEVEL_COLOUR[lev], lw=1.4, ms=5,
                label={"B1": "L9 only", "B2-lin": "L8 → L9", "B3-lin": "L7 → L8 → L9"}[v])
    ax.set_xscale("log"); ax.set_xticks([1000, 2000, 5000]); ax.set_xticklabels(["1k", "2k", "5k"]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_title(TITLE[q]); ax.set_xlabel("training cost (core-h)"); ax.set_ylim(bottom=0)
axs[0].set_ylabel("held-out L9 RMSE (ln)")
axs[-1].legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
fig.tight_layout(); fig.savefig(out / "err_vs_coreh.png", dpi=200, bbox_inches="tight")
print("L8 rows", n8, "L7 rows", n7)
for q in TITLE:
    for v in STY:
        s = a[(a.qoi == q) & (a.variant == v)].groupby("n9").rmse.median()
        print(q, v, [(int(n), round(cost(v, q, n)), round(r, 3)) for n, r in s.items()])
