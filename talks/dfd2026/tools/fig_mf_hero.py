"""Hero figure for the multi-fidelity method, on OUR data (single panel): the 95th-percentile wall shear versus rocking speed at theta = 7 deg.
 L8 (cheap): many runs, their surrogate (dashed light curve); L9: fewer runs, surrogate (dashed mid curve); L10 (expensive): few runs, chained prediction + 95 % band.
Same models and data as scripts/fig_mf_grid.py (chain L8 KRR -> L9 GP -> L10 GP on log QoI, L9 uncertainty propagated).
Usage: uv run python fig_mf_hero.py <version>   -> ../public/mf_hero_<version>.png (transparent) and a white-background copy in the scratchpad
"""
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO)]
from scripts import fig_mf_grid as g  # noqa: E402

warnings.filterwarnings("ignore")
VER = sys.argv[1] if len(sys.argv) > 1 else "v1"
OUT = Path(__file__).resolve().parent.parent / "public"
MUTE, INK, ACC = "#7a8797", "#14202e", "#e8590c"
C8, C9, C10 = "#9fc3e8", "#4f93d6", "#1d2a78"
Q, TH = "tau_95_t", 7.0

a8, a9, a10, hf = g.hydro_data(Q)
pred = g.models("chain", a8, a9, a10)
low, (mean, lo, hi), l8curve = pred(TH)
fig, ax = plt.subplots(figsize=(6.0, 4.0))
sel = lambda a: (a[0][:, 1] == TH)  # noqa: E731
ax.plot(g.GRID, l8curve, color=C8, lw=2.2, ls=(0, (1, 1.6)), zorder=1)
ax.plot(g.GRID, low, color=C9, lw=2.0, ls=(0, (4, 2)), zorder=1)
ax.fill_between(g.GRID, lo, hi, color=C10, alpha=0.14, lw=0)
ax.plot(g.GRID, mean, color=C10, lw=2.4, zorder=2)
ax.plot(a8[0][sel(a8), 0], a8[1][sel(a8)], "o", color=C8, ms=4.5, mew=0, zorder=3)
ax.plot(a9[0][sel(a9), 0], a9[1][sel(a9)], "o", color=C9, ms=8, mew=0, zorder=4)
ax.plot(a10[0][sel(a10), 0], a10[1][sel(a10)], "o", color=C10, ms=11, mew=0, zorder=5)
ax.set_yscale("log"); ax.set_xlim(14, 39)
ax.set_xlabel("rocking speed  [rpm]", color=MUTE, fontsize=13); ax.set_ylabel(r"wall shear stress $\tau_{95}$  [Pa]", color=MUTE, fontsize=13)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ax.spines[sp].set_color(MUTE)
ax.tick_params(which="both", colors=MUTE, labelcolor=MUTE, labelsize=11); ax.patch.set_alpha(0)
# direct labels, one per level, placed at the right end of each curve
xl = g.GRID[-1]
for yv, txt, col, dy in ((l8curve[-1], "L8", C8, 0), (low[-1], "L9", C9, 9), (mean[-1], "L10", C10, -9)):
    ax.annotate(txt, (xl, yv), xytext=(6, dy), textcoords="offset points", color=col, fontsize=13, va="center", fontweight="bold")
fig.tight_layout()
fig.savefig(OUT / f"mf_hero_{VER}.png", dpi=200, transparent=True)
fig.savefig(f"/oscar/scratch/eaguerov/tmp/mf_hero_{VER}_white.png", dpi=200, facecolor="white")
print("n8 n9 n10 at theta 7:", int(sel(a8).sum()), int(sel(a9).sum()), int(sel(a10).sum()))
