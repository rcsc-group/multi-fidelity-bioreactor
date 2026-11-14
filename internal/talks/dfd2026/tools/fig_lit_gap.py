"""Literature gap as a dot matrix: rows = studies, columns = what each one does. filled = yes, half = partly, hollow = no.
Every cell rests on the full-text reading of 2026-11-02 (see diary): page refs in the comments.
Usage: uv run python tools/fig_lit_gap.py -> ../public/lit_gap.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge

OUT = Path(__file__).resolve().parent.parent / "public" / "lit_gap.png"
INK, MUTE, OURS = "#14202e", "#7a8797", "#e8590c"
plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "sans-serif"})

COLS = ["rocking bag", "resolved two-phase", r"$\tau,\ k_La,\ t_{mix}$", "grid levels as data", "error bars"]
# 1 yes, 0.5 partly, 0 no
ROWS = [
    ("Zhan 2019", [1, 1, 0.5, 0, 0]),             # Fluent VOF, no turb. model; tau + exp. kLa, mixing only via vorticity; 1.7 M cells, no mesh study
    ("Seidel 2022", [1, 1, 0.5, 0, 0]),           # OpenFOAM interFoam PLIC; kLa + t_mix, no tau (Kolmogorov); GCI study, then one mesh
    ("Kim 2025", [1, 1, 1, 0, 0]),                # Basilisk VOF; tau, kLa, t_mix; convergence study, one grid used
    ("Cantarero-Rivera 2024", [0, 0, 0.5, 1, 0]),  # spinner flask, frozen rotor RANS; coarse->fine ANN correction, no UQ
    ("Kiviat 2025", [0, 0, 1, 0, 0]),             # STR 200 L-200 m3, Euler-Euler RANS k-eps; shear, kLa, t_mix
    ("Martens 2025", [0, 0, 0, 0.5, 1]),          # MF GP over culture scales (ODE simulator, no CFD); GP gives uncertainty
    ("this work", [1, 1, 1, 1, 1]),
]

fig, ax = plt.subplots(figsize=(6.4, 4.3))
R = 0.28
for i, (name, v) in enumerate(ROWS):
    y = len(ROWS) - 1 - i
    c = OURS if name == "this work" else INK
    for j, val in enumerate(v):
        ax.add_patch(plt.Circle((j, y), R, fc="white", ec=c, lw=1.6, zorder=2))
        if val == 1:
            ax.add_patch(plt.Circle((j, y), R, fc=c, ec=c, lw=1.6, zorder=3))
        elif val == 0.5:
            ax.add_patch(Wedge((j, y), R, 90, 270, fc=c, ec=c, lw=1.6, zorder=3))
    ax.text(-0.55, y, name, ha="right", va="center", fontsize=13, color=c, fontweight="bold" if name == "this work" else "normal")
for j, col in enumerate(COLS):
    ax.text(j - 0.1, len(ROWS) - 0.45, col, ha="left", va="bottom", fontsize=12, color=MUTE, rotation=35, rotation_mode="anchor")
ax.axhline(0.5, color=MUTE, lw=0.6, xmin=0.0, xmax=1.0)
ax.set_xlim(-0.5, len(COLS) - 0.5)
ax.set_ylim(-0.5, len(ROWS) + 1.6)
ax.set_aspect("equal")
ax.axis("off")
fig.savefig(OUT, dpi=220, transparent=True, bbox_inches="tight")
print("wrote", OUT)
