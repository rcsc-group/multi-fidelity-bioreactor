"""Settling-time study analysis: for each (level, out-point), compare the
checkpoint-transition's tau_mean trajectory against its INDEPENDENT
reference cold start (not a self-referential tail-stability check) --
diary.md 2026-09-06. Reports t_settle in cycles-since-checkpoint: the
first cycle after which the transition's per-cycle peak tau_mean stays
within `TOL` of the reference's own per-cycle peak for the rest of the
recorded window.

TOL=0.05 (5%) -- comfortably above the ~0.14% mean / 0.56% max rank-
induced noise floor measured on tau_mean (2026-09-06 rank-invariance
check), so a "not settled" call reflects real transient behavior, not
numerical noise.

Usage:
    uv run python scripts/analyze_settling_study.py
"""
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
MANIFEST = json.loads((ROOT / "experiments" / "settling_study_manifest.json").read_text())

BASELINE = MANIFEST["baseline"]
LEVELS = MANIFEST["levels"]
TOL = 0.05

OUT_POINTS = [
    {"rpm": 17.5, "theta": 2.0}, {"rpm": 17.5, "theta": 4.0}, {"rpm": 17.5, "theta": 7.0},
    {"rpm": 32.5, "theta": 2.0}, {"rpm": 32.5, "theta": 4.0},
    {"rpm": 37.5, "theta": 2.0}, {"rpm": 37.5, "theta": 4.0}, {"rpm": 37.5, "theta": 7.0},
]


def t_per_nd(rpm: float, theta: float) -> float:
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


def point_tag(p):
    return f"rpm{p['rpm']:g}_th{p['theta']:g}"


def t_settle(transition_peaks, reference_peaks, tol=TOL):
    n = min(len(transition_peaks), len(reference_peaks))
    reldiff = np.abs(transition_peaks[:n] - reference_peaks[:n]) / np.abs(reference_peaks[:n])
    for start in range(n):
        if np.all(reldiff[start:] < tol):
            return start, reldiff
    return None, reldiff


print(f"{'L':>2} {'rpm':>6} {'theta':>6} {'d_rpm':>7} {'d_theta':>8} {'t_settle(cyc)':>14}")
results = []
for level in LEVELS:
    lvl = str(level)
    for p in OUT_POINTS:
        tag = point_tag(p)
        transition_run_id = MANIFEST["transitions"][lvl][tag]
        reference_run_id = f"settling_ref_L{level}_{tag}"

        T_nd = t_per_nd(p["rpm"], p["theta"])

        tr = np.loadtxt(RUNS / transition_run_id / "shear_stress.dat", skiprows=1)
        ref = np.loadtxt(RUNS / reference_run_id / "shear_stress.dat", skiprows=1)

        # transition: t0 = its own t_checkpoint (first row's t minus a tiny
        # epsilon won't work -- read the actual checkpoint t from params.json)
        params = json.loads((RUNS / transition_run_id / "params.json").read_text())
        t_ckpt = params["t_checkpoint"]

        n_cycles_tr = int((tr[-1, 1] - t_ckpt) / T_nd)
        n_cycles_ref = int(ref[-1, 1] / T_nd)
        n_cycles = min(n_cycles_tr, n_cycles_ref)

        tr_peaks = per_cycle_peak(tr[:, 1], tr[:, 5], t_ckpt, T_nd, n_cycles)
        ref_peaks = per_cycle_peak(ref[:, 1], ref[:, 5], 0.0, T_nd, n_cycles)

        settle_cycle, reldiff = t_settle(tr_peaks, ref_peaks)
        d_rpm = p["rpm"] - BASELINE["rpm"]
        d_theta = p["theta"] - BASELINE["theta"]
        results.append({
            "level": level, "rpm": p["rpm"], "theta": p["theta"],
            "d_rpm": d_rpm, "d_theta": d_theta, "t_settle": settle_cycle,
            "n_cycles_available": n_cycles,
        })
        settle_str = str(settle_cycle) if settle_cycle is not None else "NEVER"
        print(f"{level:>2} {p['rpm']:>6.1f} {p['theta']:>6.1f} {d_rpm:>7.1f} {d_theta:>8.1f} {settle_str:>14}")

out_path = ROOT / "experiments" / "settling_study_results.json"
out_path.write_text(json.dumps(results, indent=2))
print(f"\nSaved {out_path}")
