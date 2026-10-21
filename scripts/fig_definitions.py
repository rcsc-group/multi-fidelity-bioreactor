"""How kLa and the mixing time are defined (one L8 run, 25 rpm, 7 deg): C*(t) with the whole-period fits that give kLa at C* = 0.10, 0.25, 0.50, and chi(t) with the crossings that give the mixing times.
C*(t) = oxygen in the liquid / saturation (scripts/postprocess.py _compute_c_star); chi(t) = 1 - sigma^2(t)/sigma^2(0) of the passive tracer in the liquid (_compute_mixing_metrics).
kLa(C* = x): slope of ln(1 - C*) vs t over one rocking period centred on the first crossing of x (_kla_period_at_threshold). dtmix(chi = y): first time chi reaches y.
Usage: uv run python scripts/fig_definitions.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402
from scripts import postprocess as pp  # noqa: E402
from scripts.test_s_estimator import curves  # noqa: E402


def main():
    run = "kmix_l8_rpm25"
    d = ROOT / "runs" / run
    p = json.loads((d / "params.json").read_text())
    ts, chi, cst = curves(run)
    T_bio, T_nd = pp._t_scales(p)
    per = T_nd * T_bio
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8))
    b = fs.level_colour(8)
    ax = axes[0]
    ax.plot(ts, cst, color=b, lw=0.9)
    for x, mk in ((0.10, fs.MK["cstar_10"]), (0.25, fs.MK["cstar_25"]), (0.50, fs.MK["cstar_50"])):
        k = int(np.argmax(cst >= x))
        m = (ts >= ts[k] - per / 2) & (ts <= ts[k] + per / 2)
        ax.plot(ts[m], cst[m], color=b, lw=3.0)
        ax.plot(ts[k], cst[k], ls="none", marker=mk, mfc="w", mec="0.15", ms=8, mew=1.3)
    ax.set_xlim(0, 350)
    ax.set_xlabel("time after release (s)")
    ax.set_ylabel(r"$C^*$ (oxygen in the liquid / saturation)")
    ax.grid(**fs.GRID_KW)
    ax = axes[1]
    ax.plot(ts, chi, color=b, lw=0.9)
    for y, mk in ((0.50, fs.MK["chi_0.50"]), (0.75, fs.MK["chi_0.75"]), (0.95, fs.MK["chi_0.95"])):
        k = int(np.argmax(chi >= y))
        ax.plot(ts[k], chi[k], ls="none", marker=mk, mfc="w", mec="0.15", ms=8, mew=1.3)
        ax.axhline(y, color=fs.GUIDE, lw=1.0, ls=":")
    ax.set_xlim(0, 350)
    ax.set_xlabel("time after release (s)")
    ax.set_ylabel(r"$\chi$ (tracer mixing, 0 = segregated, 1 = uniform)")
    ax.grid(**fs.GRID_KW)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / f"experiments/multifidelity/fig_definitions.{ext}", dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    main()
