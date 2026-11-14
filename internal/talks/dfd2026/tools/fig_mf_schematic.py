"""Multi-fidelity schematic in the grammar of Yi et al. 2407.15110 Fig. 1b (conceptual, synthetic curves, text-free; the slide sets the equation in matching colours).
 low-fidelity samples (many, light) -> low-fidelity surrogate f_L (dashed) -> linear transfer rho f_L + b (grey)
 -> high-fidelity samples (few, dark) with residual sticks (orange) = delta, the part the Gaussian process models -> f_H (solid).
Usage: uv run python fig_mf_schematic.py   -> ../public/mf_schematic.png (transparent)
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent.parent / "public"
MUTE, ACC = "#7a8797", "#e8590c"
C_L, C_H, C_T = "#7fb0e0", "#1d2a78", "#9aa5b4"

x = np.linspace(0, 1, 300)
fL = lambda u: 1.6 + 1.1 * u + 0.55 * np.sin(6.5 * u)                      # noqa: E731  low-fidelity surrogate
gT = lambda u: 0.62 * fL(u) + 0.30                                         # noqa: E731  linear transfer
fH = lambda u: gT(u) + 0.55 * np.sin(7 * u)               # noqa: E731  high-fidelity truth

xl = np.linspace(0.04, 0.96, 12)
xh = np.array([0.10, 0.30, 0.52, 0.76, 0.97])
fig, ax = plt.subplots(figsize=(5.4, 3.6))
ax.plot(x, fL(x), color=C_L, lw=2.0, ls=(0, (4, 2.5)))
ax.plot(x, gT(x), color=C_T, lw=2.2)
ax.plot(x, fH(x), color=C_H, lw=2.6)
for xi in xh:
    ax.plot([xi, xi], [gT(xi), fH(xi)], color=ACC, lw=3.0, solid_capstyle="round", zorder=3)
ax.plot(xl, fL(xl), "o", color=C_L, ms=7, mew=0, zorder=4)
ax.plot(xh, fH(xh), "o", color=C_H, ms=12, mew=0, zorder=5)
ax.set_xlim(-0.02, 1.02); ax.set_ylim(0.2, 3.7)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(MUTE); ax.spines[s].set_linewidth(1.3)
ax.set_xticks([]); ax.set_yticks([]); ax.patch.set_alpha(0)
fig.tight_layout()
fig.savefig(OUT / "mf_schematic.png", dpi=200, transparent=True)
print("wrote mf_schematic.png")
