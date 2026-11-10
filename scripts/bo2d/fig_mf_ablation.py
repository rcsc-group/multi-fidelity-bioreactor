"""Held-out L9 RMSE vs number of L9 training sites, per model variant (mf_ablation.py output).
Grammar (scripts/figstyle.py): colour = lowest level the model uses (L9 green = L9 only, L8 blue = 2 levels,
L7 sky = 3 levels); solid = upstream mfbml (ML hyperparameters); dotted = ours (fully Bayesian); circle = linear basis,
square = quadratic basis. Line = median over 20 site draws. No annotations.
  uv run python scripts/bo2d/fig_mf_ablation.py OUT.png LOG [LOG ...]"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts import figstyle as fs  # noqa: E402

rows = []
for f in sys.argv[2:]:
    rows += [l.strip().split(",") for l in open(f) if l[:3] in ("kLa", "dtm", "tau")]
d = pd.DataFrame(rows, columns="qoi,n9,rep,variant,rmse,cover90,seconds".split(","))
d[["n9", "rmse", "cover90"]] = d[["n9", "rmse", "cover90"]].astype(float)
d.to_csv(Path(__file__).parent / "state" / "mf_ablation.csv", index=False)
STYLE = {"U1": (9, "-", "o"), "B1": (9, ":", "o"),
         "U2-lin": (8, "-", "o"), "U2-quad": (8, "-", "s"), "B2-lin": (8, ":", "o"), "B2-quad": (8, ":", "s"),
         "U3-lin": (7, "-", "o"), "U3-quad": (7, "-", "s"), "B3-lin": (7, ":", "o"), "B3-quad": (7, ":", "s")}
TITLE = {"kLa_1T_10": "kLa", "dtmix_0.25": "mixing time", "tau_95": "shear stress"}
plt.rcParams.update(fs.rcparams())
fig, axs = plt.subplots(1, 3, figsize=(10, 3.2), sharey=False)
for ax, q in zip(axs, TITLE):
    g = d[d.qoi == q].groupby(["variant", "n9"]).rmse.median().unstack(0)
    for v, (lev, ls, mk) in STYLE.items():
        if v in g:
            ax.plot(g.index, g[v], ls=ls, marker=mk, color=fs.LEVEL_COLOUR[lev], mfc="none" if mk == "s" else None,
                    lw=1.4, ms=5, label=v)
    ax.set_title(TITLE[q]); ax.set_xlabel("L9 training sites"); ax.set_xticks([4, 8, 16]); ax.set_ylim(bottom=0)
axs[0].set_ylabel("held-out L9 RMSE (ln)")
axs[-1].legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=200, bbox_inches="tight")
print("median RMSE at n9 = 8:"); print(d[d.n9 == 8].groupby(["qoi", "variant"]).agg(rmse=("rmse", "median"), cover=("cover90", "median")).round(3).unstack(0).to_string())
