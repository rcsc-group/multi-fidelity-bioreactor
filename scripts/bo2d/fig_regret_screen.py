"""Regret vs core-h, v2 benchmark (2026-11-11): screen-then-L9 (synth_screen.py, N = 40 L7+L8 screen) vs L9 only,
same frozen v2 truth, seeds 21-23, kLa and t_mix side by side. Uses the logged cumulative core-h of each arm.
Grammar (scripts/figstyle.py): colour = level (screen arm drawn in L8 blue, L9-only in L9 green); thin = one seed;
thick dotted = median over seeds (only where every seed has data); x log. No annotations.
  uv run python scripts/bo2d/fig_regret_screen.py OUT.png
"""
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

LOGS = {"kla": {"screen40": [8836090, 8836093, 8836096], "l9only": [8836089, 8836092, 8836095]},
        "mix": {"screen40": [8836099, 8836102, 8836105], "l9only": [8836098, 8836101, 8836104]}}
rx = re.compile(r"^(screen40|l9only)\s+seed (\d+) round \d+ core-h\s+(\d+) regret\s+([-\d.]+)% feasible (\w+)")
TITLE = {"kla": r"$k_La$", "mix": r"$t_{\mathrm{mix}}$"}
plt.rcParams.update(fs.rcparams())
plt.rcParams.update({"font.size": 15, "axes.labelsize": 16, "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 15})
fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), sharey=True)
xs = np.geomspace(1000, 45000, 300)
for ax, obj in zip(axes, ("kla", "mix")):
    for arm, col, lab in (("screen40", fs.LEVEL_COLOUR[8], "L7 + L8 screen, then L9"), ("l9only", fs.LEVEL_COLOUR[9], "L9 only")):
        ys = []
        for j in LOGS[obj][arm]:
            pre = "bo_screen" if arm == "screen40" else "bo_real"
            pts = []
            for line in open(ROOT / "logs" / f"{pre}_{j}.out"):
                m = rx.match(line)
                if m:
                    pts.append((float(m.group(3)), max(float(m.group(4)), 0.0) if m.group(5) == "True" else 100.0))
            c, r = np.array(pts).T
            r = np.fmin.accumulate(r)
            ax.step(c, r, where="post", color=col, lw=0.8, alpha=0.35)
            i = np.searchsorted(c, xs, side="right") - 1
            ys.append(np.where(i >= 0, r[np.clip(i, 0, None)], np.nan))   # a seed stops at its HIT: carry its
            #   last (best) value forward so the median spans the full axis
        Y = np.array(ys); ok = np.isfinite(Y).all(axis=0)
        ax.step(xs[ok], np.median(Y[:, ok], axis=0), where="post", color=col, lw=2.4, ls=":", label=lab)
    ax.axhline(5, color="0.5", lw=0.8, ls="--")
    ax.set_xscale("log"); ax.set_xlabel("cost (core-hours)"); ax.set_title(TITLE[obj])
    ax.set_ylim(0, 40); ax.grid(**fs.GRID_KW)
axes[0].set_ylabel("regret (%)")
axes[1].legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=200, bbox_inches="tight")
print("wrote", sys.argv[1])
