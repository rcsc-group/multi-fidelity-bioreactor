"""Design space of a rocking bag (Tufts 'The Tradeoff' style), in (rocking speed, angle).
Boundaries are a SKETCH (artist rendition, not computed): shear damage at high speed x angle,
oxygen (kLa C >= OUR) and mixing (t_mix < 1/kLa) at low speed x angle.
Dots: conditions that earlier rocking-bag CFD studies ran (from the papers, read in full):
  Zhan 2019 Fig. 8: 15, 22, 30 rpm x 4, 5.5, 7 deg
  Seidel 2022 Fig. 10: 5, 15, 25 rpm x 5, 11.5, 16.5 deg (corners + centre)
  Kim 2025 Figs. 11-12: 15-37.5 rpm at 7 deg; 2-7 deg at 32.5 rpm
Usage: uv run python tools/fig_design_space.py -> ../public/design_space.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent.parent / "public" / "design_space.png"
INK, MUTE, WIN = "#14202e", "#7a8797", "#e8590c"
plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "sans-serif"})

fig, ax = plt.subplots(figsize=(5.6, 4.2))
x = np.linspace(3, 40, 400)
shear = 260 / x          # sketch: damage above f*theta ~ const
oxy = 70 / x             # sketch: too little oxygen below
mix = 95 / x ** 1.1 + 1.2  # sketch: mixing slower than 1/kLa below
low = np.maximum(oxy, mix)
ax.fill_between(x, low, shear, where=shear > low, color=WIN, alpha=0.13, lw=0)
for y, lab, xy in [(shear, r"$\tau = \tau_{\mathrm{damage}}$", (30, 260 / 30 + 0.6)),
                   (oxy, r"$k_La\,C_{\mathrm{sat}} = \mathrm{OUR}$", (5.2, 3.2)),
                   (mix, r"$t_{\mathrm{mix}} = 1/k_La$", (22, 95 / 22 ** 1.1 + 1.2 - 1.6))]:
    ax.plot(x, y, color=INK, lw=1.4)
    ax.text(*xy, lab, fontsize=13, color=INK)
ax.text(33, 14.5, "shear damage", color=MUTE, fontsize=12, ha="center")
ax.text(8.5, 1.0, "too little O$_2$, slow mixing", color=MUTE, fontsize=12)
ax.text(17, 8.6, "design space", color=WIN, fontsize=14, fontweight="bold", rotation=-24)

studies = {"Zhan 2019": ("s", [(r, a) for r in (15, 22, 30) for a in (4, 5.5, 7)]),
           "Seidel 2022": ("D", [(5, 5), (5, 16.5), (25, 5), (25, 16.5), (15, 11.5)]),
           "Kim 2025": ("o", [(r, 7) for r in np.arange(15, 38, 2.5)] + [(32.5, a) for a in (2, 3, 4, 5, 6)])}
for name, (m, pts) in studies.items():
    p = np.array(pts, float)
    ax.plot(p[:, 0], p[:, 1], m, ms=6, mfc="white", mec=INK, mew=1.3, ls="none", label=name)
ax.set_xlim(3, 40); ax.set_ylim(0, 18)
ax.set_xlabel("rocking speed (rpm)", color=MUTE, fontsize=13)
ax.set_ylabel(r"angle $\theta$ (deg)", color=MUTE, fontsize=13)
ax.set_xticks([]); ax.set_yticks([])
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ax.spines[sp].set_color(MUTE)
ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=12)
fig.savefig(OUT, dpi=220, transparent=True, bbox_inches="tight")
print("wrote", OUT)
