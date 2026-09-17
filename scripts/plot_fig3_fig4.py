"""Replicas of Kim et al. (2024) Figs. 3 and 4 -- the vortical structure.

Fig. 3: instantaneous vorticity overlaid with streamlines at four phases of
one rocking cycle, t/T_p = 80.00, 80.25, 80.50, 80.75.

Fig. 4: the STEADY STREAMING -- the vorticity of the time-averaged flow, and
the volume fraction, each overlaid with streamlines of the mean flow.

The distinction between the two is the whole point of Fig. 4 and is easy to
lose: Fig. 3's vorticity is curl(u) at an instant, Fig. 4's is curl(ubar) with
the average INSIDE the curl. A time average of |curl u| is a different and
much larger quantity -- the error that made an earlier settling proxy read
converged while the mean flow was still drifting 18%.

Both panels read fields the solver writes on a uniform grid:
fields/snap_*.bin for Fig. 3 and streaming_field.bin for Fig. 4,
which carries ubar itself because streamlines cannot be drawn from
streaming.dat's single spatial mean.

Fig. 3's four phases are matched on PHASE -- t mod T_p -- not on absolute
time. The recording interval is deliberately not a divisor of the period, so
frames drift across phase and a given phase is generally sampled in some other
cycle than the one Kim labels. For a settled periodic flow the two are the
same field, and matching on phase is what lets the run record five frames per
period instead of the thirteen an absolute-time match would need. The cycle
each frame actually came from is printed on its panel.

Streamlines are drawn only where there is water for essentially the whole
averaging window; elsewhere ubar is a blend of water and air and its
streamlines are an artifact of the interface sweeping the cell.

Usage:  uv run python scripts/plot_fig3_fig4.py <run_id>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
OUT_DIR = ROOT / "experiments/kimetal2024/figure_replicas"
sys.path.insert(0, str(ROOT))
from scripts.fields import (  # noqa: E402
    cell_centres, liquid_mask, load_snapshot, load_streaming_field, snapshots,
)
from scripts.postprocess import _t_scales  # noqa: E402

plt.rcParams.update({
    "mathtext.fontset": "cm", "font.family": "serif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
})

PHASES = [80.00, 80.25, 80.50, 80.75]
BAG_HALF_HEIGHT = 0.143          # b/L, geometry.b / geometry.a


def _bag_window(ax):
    """Kim crops to the bag; the domain is mostly empty above and below it."""
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(-BAG_HALF_HEIGHT, BAG_HALF_HEIGHT)
    ax.set_aspect("equal")
    ax.set_xticks([-0.5, 0.0, 0.5])
    ax.set_yticks([-0.1, 0.0, 0.1])
    ax.tick_params(which="both", direction="in", labelsize=8)


def _streamlines(ax, X, Y, ux, uy, mask, colour="k"):
    u = np.where(mask, ux, np.nan)
    v = np.where(mask, uy, np.nan)
    ax.streamplot(X[0], Y[:, 0], u, v, color=colour, linewidth=0.5,
                  density=1.1, arrowsize=0.6)


def _nearest_phase(run_dir: Path, phase: float, T_per_nd: float):
    """The recorded frame whose PHASE is closest to `phase` (mod 1)."""
    paths = snapshots(run_dir)
    if not paths:
        return None, None
    times = np.array([load_snapshot(p)["t"] for p in paths]) / T_per_nd
    d = np.abs(((times - phase + 0.5) % 1.0) - 0.5)
    k = int(np.argmin(d))
    return paths[k], float(d[k])


def fig3(run_dir: Path, T_per_nd: float) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(5.4, 7.2))
    drawn = 0
    for ax, phase in zip(axes, PHASES):
        p, err = _nearest_phase(run_dir, phase, T_per_nd)
        if p is None:
            ax.text(0.5, 0.5, "no snapshots recorded", transform=ax.transAxes,
                    ha="center", fontsize=9)
            _bag_window(ax)
            continue
        s = load_snapshot(p)
        got = s["t"] / T_per_nd
        if err > 0.03:
            print(f"  WARNING: closest recorded phase to t/T_p={phase} is off "
                  f"by {err:.3f} of a period -- cadence too coarse")
        X, Y = cell_centres(s["n"])
        m = liquid_mask(s)
        w = np.where(m, s["omega"], np.nan)
        lim = float(np.nanpercentile(np.abs(w), 99)) or 1.0
        ax.pcolormesh(X, Y, w, cmap="RdBu_r", vmin=-lim, vmax=lim,
                      shading="auto", rasterized=True)
        _streamlines(ax, X, Y, s["ux"], s["uy"], m)
        ax.set_ylabel(rf"$t/T_p={got:.2f}$", fontsize=9)   # actual cycle
        _bag_window(ax)
        drawn += 1
    axes[-1].set_xlabel(r"$x'/L_b$", fontsize=11)
    fig.tight_layout()
    out = OUT_DIR / "replicated_Fig3.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"saved {out}  ({drawn}/{len(PHASES)} phases)")


def fig4(run_dir: Path) -> None:
    path = run_dir / "streaming_field.bin"
    if not path.exists():
        print(f"  no streaming_field.bin in {run_dir.name} -- Fig 4 needs a "
              f"run recorded with the 2026-09-16 solver or later")
        return
    d = load_streaming_field(path)
    X, Y = cell_centres(d["n"])
    # Water for essentially the whole window. ubar elsewhere is a water/air
    # blend and its streamlines follow the interface, not the mean flow.
    m = d["f_acc"] > 0.99 * d["window_dt"]
    f_bar = d["f_acc"] / d["window_dt"]

    fig, axes = plt.subplots(2, 1, figsize=(5.4, 4.0))
    w = np.where(m, d["omega_bar"], np.nan)
    lim = float(np.nanpercentile(np.abs(w), 99)) or 1.0
    im = axes[0].pcolormesh(X, Y, w, cmap="RdBu_r", vmin=-lim, vmax=lim,
                            shading="auto", rasterized=True)
    fig.colorbar(im, ax=axes[0], label=r"$\bar\xi_b'$", pad=0.02)
    axes[1].pcolormesh(X, Y, f_bar, cmap="Blues", vmin=0, vmax=1,
                       shading="auto", rasterized=True)
    for ax, letter in zip(axes, "ab"):
        _streamlines(ax, X, Y, d["ubar_x"], d["ubar_y"], m)
        _bag_window(ax)
        ax.text(-0.13, 1.02, rf"$({letter})$", transform=ax.transAxes,
                fontsize=13, style="italic")
    axes[1].set_xlabel(r"$x'/L_b$", fontsize=11)
    fig.tight_layout()
    out = OUT_DIR / "replicated_Fig4.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"saved {out}  (window {d['window_dt']:.4f}, t={d['t']:.4f})")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    run_dir = ROOT / "runs" / sys.argv[1]
    params = json.loads((run_dir / "params.json").read_text())
    _, T_per_nd = _t_scales(params)
    fig3(run_dir, T_per_nd)
    fig4(run_dir)


if __name__ == "__main__":
    main()
