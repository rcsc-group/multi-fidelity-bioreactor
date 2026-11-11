"""Learning curves (state/learning_curve_<q>.json from learning_curve.py): held-out error in % (100 x RMSE ln) vs
number of runs. Top row: L7 map calibrated to L9 (a + b * map, LOO over the 33 L9 sites; % = ratio x LOO mean-only
error) vs N7. Bottom row: L9-only GP vs N9 (held-out = the other L9 sites). Line = mean over draws, dots = draws.
  uv run python scripts/bo2d/fig_learning_curve.py -> figures/bo2d/fig_learning_curve.png
"""
import json
import sys
from pathlib import Path

import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

S = ROOT / "scripts/bo2d/state"
QS = (("tau_95", r"$\tau_{95}$"), ("dtmix_0.25", r"$t_{\mathrm{mix}}$"), ("kLa_1T_10", r"$k_La$"))


def loo_mean_err(q):
    d = pd.read_csv(S / f"truth_snapshot_v2_{q}.csv"); y = np.log(d[d.level == 9].y.astype(float))
    y = y[np.isfinite(y)].to_numpy()
    return float(np.sqrt(np.mean([(y[i] - np.delete(y, i).mean()) ** 2 for i in range(len(y))])))


plt.rcParams.update(fs.rcparams())
fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.6), sharey=True)
for j, (q, lab) in enumerate(QS):
    o = json.loads((S / f"learning_curve_{q}.json").read_text()); e0 = loo_mean_err(q)
    for i, (part, key, L, scale) in enumerate((("A", "l9_calib_ratio", 7, e0), ("B", "rmse_l9", 9, 1.0))):
        r = pd.DataFrame(o[part]); r["e"] = 100 * r[key] * scale; c = fs.level_colour(L)
        ax = axes[i, j]
        ax.plot(r.n, r.e, "o", ms=3, color=c, alpha=0.4, ls="none")
        g = r.groupby("n").e.mean(); ax.plot(g.index, g.values, "-o", color=c, ms=5, lw=1.6)
        ax.set_xscale("log"); ax.grid(**fs.GRID_KW); ax.set_ylim(0, 60)
        ax.set_xticks(sorted(r.n.unique())); ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter()); ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        if i == 0:
            ax.set_title(lab)
        ax.set_xlabel("L7 runs" if i == 0 else "L9 runs")
axes[0, 0].set_ylabel("L9 error, L7 map (%)"); axes[1, 0].set_ylabel("L9 error, L9 only (%)")
fig.tight_layout(); out = ROOT / "figures/bo2d/fig_learning_curve.png"; fig.savefig(out, dpi=170); print("wrote", out)
