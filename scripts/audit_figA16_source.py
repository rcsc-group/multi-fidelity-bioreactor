"""Audit which run actually backs FigA16's L6 point: does runs/health_l6
(omega_b=3.93=~37.5rpm, geometry.b=0.071=OLD pre-fix geometry) match the
diary's claim "L6 was already at the correct RPM" (2026-08-04), or does
it contradict it? (diary.md 2026-09-04)

Usage:
    uv run python scripts/audit_figA16_source.py
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

RUNS = Path(__file__).parent.parent / "runs"


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


import json

for run_id in ["health_l6", "fig_a16_l8_rpm32p5", "fig13a_l6_rpm32.5",
               "l10_kim_seg2", "l10_kim_fig8_signed"]:
    params = json.loads((RUNS / run_id / "params.json").read_text())
    rpm = params["omega_b"] * 60 / (2 * math.pi)
    T_per_st = t_per_st(params)

    normf = np.loadtxt(RUNS / run_id / "normf.dat", skiprows=1)
    t = normf[:, 1]
    ux_rms = normf[:, 7]
    uy_rms = normf[:, 11]
    t_tp = t / T_per_st

    mask = (t_tp >= 29.0) & (t_tp <= 31.0)
    n = mask.sum()
    print(f"{run_id}: rpm={rpm:.2f} geometry.b={params['geometry']['b']} "
          f"fidelity={params['fidelity']} t_end={params['t_end']} "
          f"max_t/Tp={t_tp.max():.2f} n_in_[29,31]={n}")
    if n > 0:
        print(f"  ux_rms peak in [29,31]: {ux_rms[mask].max():.4f}  "
              f"uy_rms peak: {uy_rms[mask].max():.4f}  (Kim: ~0.80 / ~0.21)")
    else:
        print(f"  run never reached t/Tp=29 (max={t_tp.max():.2f}) -- "
              f"cannot check the [29,31] window at all")
