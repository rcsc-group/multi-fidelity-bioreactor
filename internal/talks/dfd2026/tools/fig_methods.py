"""Methods slide: three TEXT-FREE panels (light mode), transparent PNGs ../public/met_<name>.png.
 bag    : the real bag outline (embedded boundary cs = 0.5 of the L8 snapshot), half full, rocked by +-theta_max about its centre
 chain  : design points (rpm x theta) that exist at each level for kLa (same data as scripts/fig_mf_grid.py): many at L8, fewer at L9, fewest at L10
 timeline: rocking angle over the 80 settling cycles, then the tracer and oxygen release (cycle 80), then the measured window
Usage: uv run python fig_methods.py
"""
import glob
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Polygon

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO)]
from scripts import fig_mf_grid as g  # noqa: E402

RUN = REPO / "runs/rec_l8_rpm32.5_th7"
OUT = Path(__file__).resolve().parent.parent / "public"
MUTE, BLUE, DEEP, ACC = "#7a8797", "#1f6fc0", "#1d2a78", "#e8590c"
LV = {8: "#9fc3e8", 9: "#4f93d6", 10: DEEP}   # lighter to darker with fidelity


def save(fig, name):
    fig.savefig(OUT / f"met_{name}.png", dpi=200, transparent=True)
    plt.close(fig)
    print("wrote", name)


def rot(xy, deg):
    a = math.radians(deg); c, s = math.cos(a), math.sin(a)
    return xy @ np.array([[c, s], [-s, c]])


def fig_bag():
    f0 = sorted(glob.glob(str(RUN / "fields/snap_*.bin")))[0]
    cs = np.array(np.memmap(f0, dtype=np.float32, mode="r", offset=16, shape=(10, 256, 256))[4])
    r = np.where((cs > 0.5).any(axis=1))[0]; c = np.where((cs > 0.5).any(axis=0))[0]
    cs = cs[r[0]:r[-1] + 1, c[0]:c[-1] + 1]
    tmp = plt.figure().add_subplot()
    seg = max(tmp.contour(np.pad(cs, 2), [0.5]).allsegs[0], key=len)
    plt.close("all")
    seg = seg.copy()                                             # contour paths are already (x, y)
    seg -= seg.mean(axis=0); seg /= np.ptp(seg[:, 0]) * 0.5      # width -> [-1, 1]
    fig, ax = plt.subplots(figsize=(4.0, 3.0))
    for d in (-7, 7):                                            # +-theta_max ghosts
        q = rot(seg, d); ax.plot(q[:, 0], q[:, 1], color=MUTE, lw=1.0, alpha=0.55)
    outline = Polygon(seg, closed=True, fc="none", ec=DEEP, lw=2.6); ax.add_patch(outline)
    water = ax.fill_between([-1.1, 1.1], -1.0, 0.0, color=BLUE, alpha=0.25, lw=0)
    water.set_clip_path(Polygon(seg, closed=True, transform=ax.transData))
    ax.add_patch(FancyArrowPatch((-0.62, -0.62), (0.62, -0.62), connectionstyle="arc3,rad=0.45", arrowstyle="<|-|>", mutation_scale=16, color=ACC, lw=2.0))
    ytop = seg[:, 1].max()
    for x in (-0.5, 0.0, 0.5):                                   # oxygen entering at the free surface
        ax.add_patch(FancyArrowPatch((x, ytop + 0.45), (x, 0.06), arrowstyle="-|>", mutation_scale=12, color=BLUE, lw=1.6))
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.0, ytop + 0.6)
    fig.subplots_adjust(0, 0, 1, 1); save(fig, "bag")


def fig_chain():
    """Chain drawn in the grammar of Yi et al. 2407.15110 Fig. 1b: the surrogate of the lower level (dashed ghost) is mapped linearly onto the next level
    (grey curve); the residual (orange sticks) is what the Gaussian process models, and the level's own surrogate (solid colour) is transfer + residual.
    Conceptual (synthetic curves), not data."""
    rng = np.random.default_rng(3)
    x = np.linspace(0, 1, 200)
    f8 = lambda u: 0.55 + 0.9 * u + 0.28 * np.sin(7 * u)                     # noqa: E731
    g9 = lambda u: 0.70 * f8(u) + 0.55                                       # noqa: E731  linear transfer of L8
    f9 = lambda u: g9(u) + 0.22 * np.sin(5 * u + 0.5)                        # noqa: E731  residual structure
    g10 = lambda u: 0.92 * f9(u) - 0.12                                      # noqa: E731
    f10 = lambda u: g10(u) + 0.20 * np.sin(6 * u + 2.0)                     # noqa: E731
    fig, axs = plt.subplots(1, 3, figsize=(9.4, 2.0), gridspec_kw=dict(wspace=0.55))
    layers = ((8, 34, None, None, f8), (9, 6, f8, g9, f9), (10, 4, f9, g10, f10))
    for ax, (lvl, n, prev, g, cur) in zip(axs, layers):
        xs = np.sort(rng.uniform(0.05, 0.95, n)) if lvl == 8 else np.linspace(0.12, 0.88, n); col = LV[lvl]
        if prev is None:
            ys = cur(xs) + rng.normal(0, 0.05, n)
        else:
            ax.plot(x, prev(x), color=LV[lvl - 1], lw=1.6, ls=(0, (3, 2)))      # surrogate of the lower level (ghost)
            ax.plot(x, g(x), color=MUTE, lw=1.8)                                # its linear transfer
            ys = cur(xs)
            for xi, yi in zip(xs, ys):
                ax.plot([xi, xi], [g(xi), yi], color=ACC, lw=2.4, solid_capstyle="round", zorder=3)   # residual
        ax.plot(x, cur(x), color=col, lw=2.2, alpha=0.9)
        ax.plot(xs, ys, "o", color=col, ms=[3.6, 7.5, 10.0][lvl - 8], mew=0, zorder=4)
        ax.set_xlim(0, 1); ax.set_ylim(0.0, 2.6)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(MUTE); ax.spines[sp].set_linewidth(1.2)
        ax.set_xticks([]); ax.set_yticks([]); ax.patch.set_alpha(0)
    for i in (0, 1):
        x0 = axs[i].get_position().x1 + 0.008; x1 = axs[i + 1].get_position().x0 - 0.008
        fig.add_artist(FancyArrowPatch((x0, 0.5), (x1, 0.5), transform=fig.transFigure, arrowstyle="-|>", mutation_scale=16, color=ACC, lw=2.2))
    fig.subplots_adjust(left=0.02, right=0.99, top=0.97, bottom=0.08)
    save(fig, "chain")


def fig_timeline():
    fig, ax = plt.subplots(figsize=(4.6, 2.2))
    t1 = np.linspace(0, 5, 900); t2 = np.linspace(7.4, 12.4, 900)
    ax.plot(t1, np.sin(2 * np.pi * t1), color=MUTE, lw=1.6)                    # settling cycles (80 in the runs; drawn broken)
    ax.plot([5.3, 6.1], [0.0, 0.0], color=MUTE, lw=1.6, ls=(0, (1, 2)))
    ax.plot(t2, np.sin(2 * np.pi * t2), color=BLUE, lw=1.8)
    ax.plot([6.75, 6.75], [-4.0, 1.5], color=ACC, lw=2.8)                      # release at cycle 80
    s = np.linspace(0, 5.6, 100)
    ax.plot(7.0 + s, -2.1 - 1.4 * (1 - np.exp(-s / 4.2)), color=BLUE, lw=2.4)               # oxygen
    ax.plot(7.0 + s, -2.1 - 1.4 * (1 - np.exp(-s / 1.8)), color=DEEP, lw=2.4, ls=(0, (4, 2)))  # tracer
    ax.set_ylim(-4.2, 1.7); ax.set_xlim(-0.2, 12.9); ax.axis("off")
    fig.subplots_adjust(0, 0, 1, 1); save(fig, "timeline")


if __name__ == "__main__":
    fig_bag(); fig_chain(); fig_timeline()
