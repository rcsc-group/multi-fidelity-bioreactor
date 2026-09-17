"""Replica of Kim et al. (2024) Fig. A.18 -- 2D against 3D.

Time evolution of the spatially averaged (a) and rms (b) liquid velocity
components in the non-inertial frame, over two cycles, for a 2D and a 3D
simulation of the same condition. Kim runs this to justify studying the
bioreactor in 2D; reproducing the justification is what lets our own 2D
results stand on the same footing.

No solver changes were needed for the 3D half. The embedded bag is a function
of x and y only, so it extrudes into a slab under -grid=octree and the liquid
volume comes out identical at both dimensionalities (0.286, measured). What
WAS missing was the spanwise component: normf.dat logged x and y and nothing
else, so the one quantity the comparison turns on was the one quantity absent.
u_z now has its own columns, written as zeros in 2D so the file keeps the same
shape at both dimensionalities.

The 3D result to look for is not agreement in u_z -- that is zero by symmetry
in the mean either way -- but whether u_z,rms stays small compared with the
in-plane components. If it does, the flow is two-dimensional in the sense that
matters and a 2D study is sound; if it does not, every 2D number in this
project inherits the error.

Usage:  uv run python scripts/plot_figA18.py <run_2d> <run_3d>
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
OUT = ROOT / "experiments/kimetal2024/figure_replicas/replicated_FigA18.png"
sys.path.insert(0, str(ROOT))
from scripts.postprocess import _t_scales  # noqa: E402

plt.rcParams.update({
    "mathtext.fontset": "cm", "font.family": "serif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
})

_COL_T = 1
_COL_RMS = {"x": 7, "y": 11, "z": 17}
_COL_SAVG = {"x": 14, "y": 15, "z": 16}
COMPONENT_COLOUR = {"x": "#0072B2", "y": "#D55E00", "z": "#009E73"}
WINDOW_CYCLES = 2.0


def load(run_id: str):
    run_dir = ROOT / "runs" / run_id
    data = np.loadtxt(run_dir / "normf.dat", skiprows=1)
    if data.shape[1] <= max(_COL_RMS.values()):
        raise SystemExit(
            f"{run_id}: normf.dat has {data.shape[1]} columns -- this run "
            f"predates the u_z columns and cannot feed Fig A.18")
    params = json.loads((run_dir / "params.json").read_text())
    _, T_per_nd = _t_scales(params)
    t_tp = data[:, _COL_T] / T_per_nd
    # The last whole two cycles, so both runs are compared over settled flow
    # rather than over whatever each happened to end on.
    hi = np.floor(t_tp.max())
    m = (t_tp >= hi - WINDOW_CYCLES) & (t_tp <= hi)
    return t_tp[m] - (hi - WINDOW_CYCLES), data[m]


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    runs = {"2D": load(sys.argv[1]), "3D": load(sys.argv[2])}

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6))
    for label, (t, d) in runs.items():
        ls = "-" if label == "2D" else "--"
        for comp in ("x", "y", "z"):
            axes[0].plot(t, d[:, _COL_SAVG[comp]], color=COMPONENT_COLOUR[comp],
                         ls=ls, lw=1.2)
            axes[1].plot(t, d[:, _COL_RMS[comp]], color=COMPONENT_COLOUR[comp],
                         ls=ls, lw=1.2)
        print(f"{label}: u_z,rms / u_x,rms = "
              f"{d[:, _COL_RMS['z']].max() / max(d[:, _COL_RMS['x']].max(), 1e-12):.4f}")

    from matplotlib.lines import Line2D
    comps = [Line2D([], [], color=c, lw=1.2, label=rf"$u'_{k}$")
             for k, c in COMPONENT_COLOUR.items()]
    dims = [Line2D([], [], color="0.3", lw=1.2, ls=ls, label=lab)
            for lab, ls in (("2D", "-"), ("3D", "--"))]
    axes[0].set_ylabel(r"$\langle u'_i\rangle/U_b$", fontsize=12)
    axes[1].set_ylabel(r"$u'_{i,rms}/U_b$", fontsize=12)
    for ax, letter in zip(axes, "ab"):
        ax.set_xlabel(r"$t/T_p$", fontsize=12)
        ax.set_xlim(0, WINDOW_CYCLES)
        ax.tick_params(which="both", direction="in")
        ax.text(-0.14, 1.03, rf"$({letter})$", transform=ax.transAxes,
                fontsize=13, style="italic")
    first = axes[1].legend(handles=comps, fontsize=8.5, frameon=False,
                           loc="upper left", bbox_to_anchor=(1.02, 1.0))
    axes[1].add_artist(first)
    axes[1].legend(handles=dims, fontsize=8.5, frameon=False,
                   loc="lower left", bbox_to_anchor=(1.02, 0.0))
    fig.tight_layout()
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
