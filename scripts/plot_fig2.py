"""Replica of Kim et al. (2024) Fig. 2 -- the rocking profile.

The rocking angle together with the spatially averaged x' and y' velocities in
the water, normalised by U_b, over two windows: the early stage of oscillation
and the cycles around t/T_p = 80, where the tracer and oxygen are released.

Both velocity series come straight from normf.dat; the angle is not logged
anywhere and does not need to be, since it is the forcing:

    theta(t) = theta_max * sin(omega_b t)

with the solver's 3-cycle soft-start ramp applied, which is why the early
window shows the amplitude growing rather than starting at its full value.
Reading the angle off a field snapshot instead would tie this figure to a
recording cadence for a quantity that is known exactly.

The velocities are the SIGNED means, normf.dat columns ux_liq_savg /
uy_liq_savg. Basilisk's normf(v).avg -- which is what ux_liq_avg is -- sums
|v|, so it rectifies the signal: it oscillates at twice the rocking frequency
and never goes negative, and the phase relationship this figure exists to show
disappears. Runs recorded before those columns existed fall back to the
rectified ones and are labelled as such, because a reader cannot tell the two
apart by eye on a single panel.

The means are over the embedded volume in units where U_b is 1, which is
already Kim's normalisation.

Usage:  uv run python scripts/plot_fig2.py [run_id]
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
OUT = ROOT / "experiments/kimetal2024/figure_replicas/replicated_Fig2.png"
sys.path.insert(0, str(ROOT))
from scripts.postprocess import _t_scales  # noqa: E402
from scripts import figstyle as fs  # noqa: E402

plt.rcParams.update({
    "mathtext.fontset": "cm", "font.family": "serif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
})

DEFAULT_RUN = "fig9_l9_rpm32.5"
N_RAMP_CYCLES = 3        # BioReactor.c: soft start over three cycles
WINDOWS = [(0.0, 6.0), (78.0, 84.0)]

_COL_T = 1
_COL_UX_AVG, _COL_UY_AVG = 6, 10        # |.| means -- every run has these
_COL_UX_SAVG, _COL_UY_SAVG = 14, 15     # signed means -- runs after 2026-09-16


def theta_deg(t_tp, theta_max: float) -> np.ndarray:
    """Rocking angle in degrees against t/T_p, including the soft start."""
    ramp = np.clip(t_tp / N_RAMP_CYCLES, 0.0, 1.0)
    return theta_max * ramp * np.sin(2 * math.pi * t_tp)


def main() -> None:
    run = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_RUN
    run_dir = ROOT / "runs" / run
    params = json.loads((run_dir / "params.json").read_text())
    T_bio, T_per_nd = _t_scales(params)
    data = np.loadtxt(run_dir / "normf.dat", skiprows=1)
    t_tp = data[:, _COL_T] / T_per_nd
    signed = data.shape[1] > _COL_UY_SAVG
    if signed:
        ux, uy = data[:, _COL_UX_SAVG], data[:, _COL_UY_SAVG]
    else:
        ux, uy = data[:, _COL_UX_AVG], data[:, _COL_UY_AVG]
    th_max = params["theta_max"][0]

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6))
    for ax, (lo, hi) in zip(axes, WINDOWS):
        m = (t_tp >= lo) & (t_tp <= hi)
        if not m.any():
            ax.text(0.5, 0.5, f"no data in $t/T_p\\in[{lo:g},{hi:g}]$",
                    transform=ax.transAxes, ha="center", fontsize=9)
            ax.set_xlim(lo, hi)
            continue
        ax.plot(t_tp[m], theta_deg(t_tp[m], th_max) / th_max, color="0.55",
                lw=1.6, label=r"$\theta_b/\theta_{b,max}$")
        bars = "" if signed else r"|"
        ax.plot(t_tp[m], ux[m], color=fs.COMPONENT_COLOUR["x"], lw=1.2,
                label=rf"$\langle {bars}u'_x{bars}\rangle/U_b$")
        ax.plot(t_tp[m], uy[m], color=fs.COMPONENT_COLOUR["y"], lw=1.2,
                label=rf"$\langle {bars}u'_y{bars}\rangle/U_b$")
        ax.set_xlim(lo, hi)
        ax.set_xlabel(r"$t/T_p$", fontsize=12)
        ax.tick_params(which="both", direction="in")
        ax.axhline(0.0, color="0.85", lw=0.8, zorder=0)
    axes[0].set_ylabel("normalized", fontsize=12)
    axes[1].legend(fontsize=9, frameon=False, loc="center left",
                   bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"{run}: {len(t_tp)} samples, t/T_p up to {t_tp.max():.1f}, "
          f"{'signed' if signed else 'RECTIFIED (pre-2026-09-16 run)'} means")
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
