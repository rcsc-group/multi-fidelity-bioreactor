"""Regenerate docs/kimetal2024/figure_replicas/replicated_FigA16_{a,b}.png
(grid convergence, u'_x,rms and u'_y,rms vs t/T_p) from a genuinely
correct L6 point (diary.md 2026-09-04).

The committed figure's L6 point came from runs/health_l6 -- despite the
2026-08-04 diary entry's claim that "L6 was already at the correct RPM,"
directly checking its params.json shows omega_b=3.93 (~37.5rpm, not
Kim's 32.5rpm baseline) AND geometry.b=0.071 (the pre-2026-08-03-fix
value). Its ux_rms peak in t/Tp=[29,31] is 0.39 -- exactly the "half of
Kim's value" signature the OLD geometry/H_bio bug produced (see
tests/verification/test_kim_fig_a16_velocity_rms.py's own docstring).
So the committed figure was comparing two DIFFERENT physical conditions,
not a genuine L6-vs-L8 grid-convergence check, and was never caught
because nothing checked the L6 point against an external reference.

Replacement L6 source: runs/fig13a_l6_rpm32.5 (this session's L6 sweep,
2026-09-01/02) -- correct RPM, correct geometry, reaches t/Tp=33 (covers
the [29,31] window). Gives ux_rms=0.774, uy_rms=0.212, matching both
Kim (~0.80/~0.21) and the L8 point (0.771/0.210) closely -- a genuine
grid-convergence result, unlike the old figure.

Usage:
    uv run python scripts/plot_figA16_current.py
"""
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RUNS = Path(__file__).parent.parent / "runs"
OUT_DIR = Path(__file__).parent.parent / "docs/kimetal2024/figure_replicas"

L6_RUN = "fig13a_l6_rpm32.5"
L8_RUN = "fig_a16_l8_rpm32p5"


def t_per_st(params):
    omega_b = params["omega_b"]
    L_bio = params["geometry"]["a"]
    H_bio = 2 * params["geometry"]["b"]
    th_max = math.radians(params["theta_max"][0])
    T_per = 2 * math.pi / omega_b
    V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(th_max))
    U_bio = V_bio / (H_bio * 0.5) / T_per
    T_bio = L_bio / U_bio
    return T_per / T_bio


def load(run_id):
    params = json.loads((RUNS / run_id / "params.json").read_text())
    normf = np.loadtxt(RUNS / run_id / "normf.dat", skiprows=1)
    t_tp = normf[:, 1] / t_per_st(params)
    mask = (t_tp >= 29.0) & (t_tp <= 31.0)
    return t_tp[mask], normf[mask, 7], normf[mask, 11]  # t/Tp, ux_rms, uy_rms


t6, ux6, uy6 = load(L6_RUN)
t8, ux8, uy8 = load(L8_RUN)

for comp, u6, u8, letter in [("x", ux6, ux8, "a"), ("y", uy6, uy8, "b")]:
    fig, ax = plt.subplots(figsize=(5.0, 4.5))
    ax.plot(t6, u6, color="crimson", lw=2.2, ls="--", label=r"$n_L=2^6$")
    ax.plot(t8, u8, color="mediumorchid", lw=1.6, ls="--", label=r"$n_L=2^8$")
    ax.set_xlabel(r"$t/T_p$", fontsize=13)
    ax.set_ylabel(rf"$\langle u'_{{{comp},rms}}\rangle/U_b$", fontsize=13)
    ax.set_ylim(0.0, 1.2)
    ax.legend(fontsize=11, loc="upper right")
    ax.text(-0.18, 1.02, rf"$({letter})$", transform=ax.transAxes, fontsize=15, style="italic")
    fig.tight_layout()
    out_path = OUT_DIR / f"replicated_FigA16_{letter}.png"
    fig.savefig(out_path, dpi=150)
    print(f"saved {out_path}")
