"""Replica of Kim et al. (2024) Fig. 1 -- the two frames of reference.

The inertial frame O, in which the bag rocks, and the non-inertial frame O' in
which it is stationary and gravity swings instead. This is the geometric
statement the whole solver rests on: BioReactor.c never moves the embedded
boundary, it rotates the body force, and Kim's Fig. 1 is the picture of why
those are the same problem.

A schematic, not a result -- every number in it is a parameter, so it is drawn
from the same params.json the runs use rather than from any output. If the
geometry in that file ever stops matching this drawing, the drawing is the one
that is wrong.

Usage:  uv run python scripts/plot_fig1.py [theta_deg]
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, FancyArrow, Rectangle

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
OUT = ROOT / "experiments/kimetal2024/figure_replicas/replicated_Fig1.png"

plt.rcParams.update({
    "mathtext.fontset": "cm", "font.family": "serif", "axes.linewidth": 1.2,
})

L_X, L_Y = 1.0, 0.286      # bag width and height, normalised by L_b
FILL = 0.5


def _bag(ax, theta_deg: float, edge: str, ls: str, label: str,
         water_level: bool = False):
    """The bag outline and its water, rotated by theta about the origin.

    `water_level` draws the surface HORIZONTAL and clipped to the tilted bag,
    which is what the inertial frame shows: there the bag rocks and gravity
    does not, so the water stays level. Rotating the water with the bag would
    draw the non-inertial picture inside the inertial panel and contradict the
    one thing the figure exists to say.
    """
    th = math.radians(theta_deg)
    R = np.array([[math.cos(th), -math.sin(th)],
                  [math.sin(th), math.cos(th)]])
    corners = np.array([[-L_X / 2, -L_Y / 2], [L_X / 2, -L_Y / 2],
                        [L_X / 2, L_Y / 2], [-L_X / 2, L_Y / 2],
                        [-L_X / 2, -L_Y / 2]])
    water = np.array([[-L_X / 2, -L_Y / 2], [L_X / 2, -L_Y / 2],
                      [L_X / 2, -L_Y / 2 + FILL * L_Y],
                      [-L_X / 2, -L_Y / 2 + FILL * L_Y],
                      [-L_X / 2, -L_Y / 2]])
    c = corners @ R.T
    ax.plot(c[:, 0], c[:, 1], color=edge, ls=ls, lw=1.6, label=label)
    if water_level:
        # Level surface at the fill height, clipped to the tilted bag.
        y_s = -L_Y / 2 + FILL * L_Y
        patch = plt.Polygon(c[:-1], closed=True, facecolor="none",
                            edgecolor="none")
        ax.add_patch(patch)
        band = ax.fill_between([-1.0, 1.0], -1.0, y_s, color="#0072B2",
                               alpha=0.18, lw=0)
        band.set_clip_path(patch)
    else:
        w = water @ R.T
        ax.fill(w[:, 0], w[:, 1], color="#0072B2", alpha=0.18, lw=0)
    return R


def _axes_pair(ax, R, colour, primed: bool):
    o = np.zeros(2)
    for vec, name in ((np.array([0.30, 0.0]), "x"), (np.array([0.0, 0.16]), "y")):
        v = R @ vec
        ax.add_patch(FancyArrow(o[0], o[1], v[0], v[1], width=0.002,
                                head_width=0.02, head_length=0.025,
                                color=colour, length_includes_head=True))
        tag = rf"${name}'$" if primed else rf"${name}$"
        ax.annotate(tag, xy=(v[0] * 1.18, v[1] * 1.18), color=colour,
                    fontsize=12, ha="center", va="center")


def panel(ax, theta_deg: float, non_inertial: bool):
    if non_inertial:
        # The bag is stationary; gravity swings by -theta instead.
        R = _bag(ax, 0.0, "#0072B2", "--", r"$O'$ (non-inertial)")
        _axes_pair(ax, np.eye(2), "#0072B2", primed=True)
        g = np.array([math.sin(math.radians(theta_deg)),
                      -math.cos(math.radians(theta_deg))]) * 0.20
        ax.add_patch(FancyArrow(0.34, 0.12, g[0], g[1], width=0.003,
                                head_width=0.022, head_length=0.028,
                                color="0.3", length_includes_head=True))
        ax.annotate(r"$\mathbf{g}'$", xy=(0.34 + g[0], 0.12 + g[1] - 0.035),
                    color="0.3", fontsize=12, ha="center")
        ax.set_title(r"$O'$: bag fixed, gravity rotates", fontsize=10)
    else:
        R = _bag(ax, theta_deg, "0.35", "-", r"$O$ (inertial)",
                 water_level=True)
        _axes_pair(ax, np.eye(2), "0.35", primed=False)
        ax.add_patch(FancyArrow(0.34, 0.12, 0.0, -0.20, width=0.003,
                                head_width=0.022, head_length=0.028,
                                color="0.3", length_includes_head=True))
        ax.annotate(r"$\mathbf{g}$", xy=(0.34, -0.12), color="0.3",
                    fontsize=12, ha="center")
        ax.add_patch(Arc((0, 0), 0.52, 0.52, theta1=0.0, theta2=theta_deg,
                         color="#D55E00", lw=1.4))
        ax.annotate(rf"$\theta_b={theta_deg:g}^\circ$", xy=(0.30, 0.045),
                    color="#D55E00", fontsize=11)
        ax.set_title(r"$O$: bag rocks, gravity fixed", fontsize=10)

    # L_x measures the BAG, so in the inertial panel it tilts with it.
    ends = np.array([[-L_X / 2, -0.24], [L_X / 2, -0.24]]) @ R.T
    ax.annotate("", xy=tuple(ends[0]), xytext=tuple(ends[1]),
                arrowprops=dict(arrowstyle="<->", color="0.45", lw=1.0))
    mid = (ends[0] + ends[1]) / 2
    ax.annotate(r"$L_x$", xy=(mid[0], mid[1] - 0.035), color="0.45",
                fontsize=11, ha="center")
    side = np.array([[0.60, -L_Y / 2], [0.60, L_Y / 2]]) @ R.T
    ax.annotate("", xy=tuple(side[0]), xytext=tuple(side[1]),
                arrowprops=dict(arrowstyle="<->", color="0.45", lw=1.0))
    smid = (side[0] + side[1]) / 2
    ax.annotate(r"$L_y$", xy=(smid[0] + 0.045, smid[1]), color="0.45",
                fontsize=11, va="center")
    ax.set_xlim(-0.72, 0.72)
    ax.set_ylim(-0.34, 0.34)
    ax.set_aspect("equal")
    ax.set_axis_off()


def main() -> None:
    theta = float(sys.argv[1]) if len(sys.argv) > 1 else 7.0
    # Drawn larger than 7 deg would read: at the real angle the two panels look
    # identical on the page, which defeats the figure. Stated, not hidden.
    theta_drawn = max(theta, 18.0)
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.4))
    panel(axes[0], theta_drawn, non_inertial=False)
    panel(axes[1], theta_drawn, non_inertial=True)
    fig.tight_layout()
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"drawn at {theta_drawn:g} deg for legibility "
          f"(the runs use {theta:g} deg)")
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
