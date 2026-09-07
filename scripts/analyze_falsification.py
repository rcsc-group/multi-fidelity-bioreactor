"""Falsification-test analysis (diary.md 2026-09-07): does the cold-start-
beats-warm-start anomaly reverse direction when warm-starting from a WEAKER
forcing (theta=2, less kinetic energy) into a STRONGER one (theta=4 or 7)?

Same rigorous methodology as compare_coldstart_vs_warmstart_settling.py:
per-cycle peak tau_mean, transition compared against an independent
reference cold start's own converged (last-5-cycle-mean) value.

Usage:
    uv run python scripts/analyze_falsification.py
"""
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
TOL = 0.05
MANIFEST = json.loads((ROOT / "experiments" / "settling_study_falsification_manifest.json").read_text())


def t_per_nd(rpm, theta):
    omega_b = rpm * 2 * math.pi / 60.0
    L_bio, H_bio = 0.25, 2 * 0.03575
    th = math.radians(theta)
    T_per = 2 * math.pi / omega_b
    V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(th))
    U_bio = V_bio / (H_bio * 0.5) / T_per
    T_bio = L_bio / U_bio
    return T_per / T_bio


def per_cycle_peak(t, tau_mean, t0, T_nd, n_cycles):
    cycles = (t - t0) / T_nd
    out = np.full(n_cycles, np.nan)
    for c in range(n_cycles):
        mask = (cycles >= c) & (cycles < c + 1)
        if mask.sum() > 0:
            out[c] = tau_mean[mask].max()
    return out


def cold_ref_run_id(level, rpm, theta):
    if (rpm, theta) == (32.5, 7.0):
        return f"settling_baseline_L{level}"
    return f"settling_ref_L{level}_rpm{rpm:g}_th{theta:g}"


def cold_start_settle(level, rpm, theta, tol=TOL, n_tail=5):
    run_id = cold_ref_run_id(level, rpm, theta)
    d = np.loadtxt(RUNS / run_id / "shear_stress.dat", skiprows=1)
    T_nd = t_per_nd(rpm, theta)
    n_cycles = int(d[-1, 1] / T_nd)
    peaks = per_cycle_peak(d[:, 1], d[:, 5], 0.0, T_nd, n_cycles)
    target = float(np.nanmean(peaks[-n_tail:]))
    reldiff = np.abs(peaks - target) / target
    for start in range(len(reldiff)):
        if np.all(reldiff[start:] < tol):
            return start
    return None


def warmstart_settle(level, to_rpm, to_theta, run_id, tol=TOL, n_tail=5):
    ref_run_id = cold_ref_run_id(level, to_rpm, to_theta)
    params = json.loads((RUNS / run_id / "params.json").read_text())
    t_ckpt = params["t_checkpoint"]
    T_nd = t_per_nd(to_rpm, to_theta)
    tr = np.loadtxt(RUNS / run_id / "shear_stress.dat", skiprows=1)
    ref = np.loadtxt(RUNS / ref_run_id / "shear_stress.dat", skiprows=1)
    n_cycles_tr = int((tr[-1, 1] - t_ckpt) / T_nd)
    n_cycles_ref = int(ref[-1, 1] / T_nd)
    tr_peaks = per_cycle_peak(tr[:, 1], tr[:, 5], t_ckpt, T_nd, n_cycles_tr)
    ref_peaks = per_cycle_peak(ref[:, 1], ref[:, 5], 0.0, T_nd, n_cycles_ref)
    target = float(np.nanmean(ref_peaks[-n_tail:]))
    reldiff = np.abs(tr_peaks - target) / target
    for start in range(len(reldiff)):
        if np.all(reldiff[start:] < tol):
            return start
    return None


print(f"{'L':>2} {'from':>10} {'to':>10} {'coldstart':>10} {'warmstart':>10} {'faster?':>8}")
for level_str, edges in MANIFEST["edges"].items():
    level = int(level_str)
    for tag, e in edges.items():
        frm, to, run_id = e["from"], e["to"], e["run_id"]
        try:
            cold = cold_start_settle(level, to["rpm"], to["theta"])
            warm = warmstart_settle(level, to["rpm"], to["theta"], run_id)
        except (FileNotFoundError, OSError):
            frm_str = f"{frm['rpm']:g}/{frm['theta']:g}"
            to_str = f"{to['rpm']:g}/{to['theta']:g}"
            print(f"{level:>2} {frm_str:>10} {to_str:>10} {'--':>10} {'--':>10} {'pending':>8}")
            continue
        cold_str = str(cold) if cold is not None else "NEVER"
        warm_str = str(warm) if warm is not None else "NEVER"
        if cold is not None and warm is not None:
            faster = "YES" if warm < cold else ("SAME" if warm == cold else "NO")
        else:
            faster = "?"
        frm_str = f"{frm['rpm']:g}/{frm['theta']:g}"
        to_str = f"{to['rpm']:g}/{to['theta']:g}"
        print(f"{level:>2} {frm_str:>10} {to_str:>10} {cold_str:>10} {warm_str:>10} {faster:>8}")
