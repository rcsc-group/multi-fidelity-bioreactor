"""Regret vs core-h for the two BO arms (2026-11-08). Synthetic held-out logs (re-costed with measured costs:
L7 5, L8 26, L9 376 incl. prerequisites; l9only 371) or, later, the real arms.
Grammar (scripts/figstyle.py): colour = level (chain drawn in L7 sky, L9-only in L9 green); thin translucent = one
seed; thick = median over seeds; dotted = ours. No annotations; claims go in the speaker notes.
  uv run python scripts/bo2d/fig_regret.py OUT.png LOG [LOG ...]"""
import re
import sys
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts import figstyle as fs  # noqa: E402

CC = {7: 5.0, 8: 26.0, 9: 376.0}; CL = 371.0
import os
OFF = float(os.environ.get("CHAIN_OFFSET", "0"))   # fair-seed charge for the chain (214 core-h, 2026-11-08)
rx = re.compile(r"^(chain|l9only)\s+seed (\d+) round \d+ core-h\s+\d+ regret\s+([-\d.]+)% feasible (\w+) .*picks L7/L8/L9 (\d+)/(\d+)/(\d+)")
cur = defaultdict(list); spent = defaultdict(float)
for f in sys.argv[2:]:
    for line in open(f):
        m = rx.match(line)
        if not m:
            continue
        arm, seed, reg, feas = m.group(1), int(m.group(2)), float(m.group(3)), m.group(4) == "True"
        n7, n8, n9 = (int(m.group(i)) for i in (5, 6, 7))
        if (arm, seed) not in spent and arm == "chain":
            spent[(arm, seed)] = OFF
        spent[(arm, seed)] += n7 * CC[7] + n8 * CC[8] + n9 * CC[9] if arm == "chain" else n9 * CL
        cur[(arm, seed)].append((spent[(arm, seed)], max(reg, 0.0) if feas else np.nan))
plt.rcParams.update(fs.rcparams())
fig, ax = plt.subplots(figsize=(5.2, 3.6))
xs = np.geomspace(100, 15000, 200)
for arm, col, lab in (("chain", fs.LEVEL_COLOUR[7], "L7 → L8 → L9"), ("l9only", fs.LEVEL_COLOUR[9], "L9 only")):
    ys = []
    for (a, s), pts in cur.items():
        if a != arm:
            continue
        c, r = np.array(pts).T
        r = np.fmin.accumulate(np.nan_to_num(r, nan=100.0))          # best recommendation so far
        ax.step(c, r, where="post", color=col, lw=0.8, alpha=0.35)
        i = np.searchsorted(c, xs, side="right") - 1                  # step: last value at or before x
        ys.append(np.where(i >= 0, r[np.clip(i, 0, None)], np.nan))
    Y = np.array(ys); ok = np.isfinite(Y).all(axis=0)                   # median only where every seed has data
    ax.step(xs[ok], np.median(Y[:, ok], axis=0), where="post", color=col, lw=2.4, ls=":", label=lab)
ax.set_xscale("log"); ax.set_xlabel("cost (core-hours)"); ax.set_ylabel("regret (%)")
ax.set_ylim(0, 60); ax.grid(**fs.GRID_KW)
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=200, bbox_inches="tight")
print("wrote", sys.argv[1])
