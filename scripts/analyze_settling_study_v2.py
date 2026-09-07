"""Settling-time analysis for the corrected nearest-neighbor-chain
transitions (diary.md 2026-09-06). Same methodology as
analyze_settling_study.py (compare tau_mean per-cycle peak against an
independent reference, not self-referential), applied to the v2 chain
topology and its 60-cycle recording window.

Usage:
    uv run python scripts/analyze_settling_study_v2.py
"""
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
MANIFEST = json.loads((ROOT / "experiments" / "settling_study_v2_manifest.json").read_text())
TOL = 0.05


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


def reference_steady_value(reference_peaks, n_tail=5):
    """Reference's own steady-state target: mean of its last n_tail cycles'
    peaks. Verified (2026-09-06) that references are stable to <1.5%
    std/mean over their own last 5 cycles, so this is a single, reusable
    target -- not requiring the transition and reference to share a cycle
    count. [FIX 2026-09-06: the original per-cycle-index comparison
    truncated to min(transition_cycles, reference_cycles), which silently
    cut off comparisons at whatever length the (shorter) reference run
    happened to be -- read as "never settles" when it was really "ran out
    of reference data to check," not a real transient.]"""
    return float(np.nanmean(reference_peaks[-n_tail:]))


def find_settle(transition_peaks, reference_peaks, tol=TOL):
    target = reference_steady_value(reference_peaks)
    reldiff = np.abs(transition_peaks - target) / abs(target)
    n = len(reldiff)
    for start in range(n):
        if np.all(reldiff[start:] < tol):
            return start, reldiff
    return None, reldiff


print(f"{'L':>2} {'from':>14} {'to':>14} {'d_rpm':>7} {'d_theta':>8} {'t_settle(cyc)':>14}")
results = []
for level_str, edges in MANIFEST["edges"].items():
    level = int(level_str)
    for tag, e in edges.items():
        frm, to, run_id = e["from"], e["to"], e["run_id"]
        ref_run_id = f"settling_ref_L{level}_rpm{to['rpm']:g}_th{to['theta']:g}"

        params = json.loads((RUNS / run_id / "params.json").read_text())
        t_ckpt = params["t_checkpoint"]

        T_nd = t_per_nd(to["rpm"], to["theta"])

        tr = np.loadtxt(RUNS / run_id / "shear_stress.dat", skiprows=1)
        ref = np.loadtxt(RUNS / ref_run_id / "shear_stress.dat", skiprows=1)

        n_cycles_tr = int((tr[-1, 1] - t_ckpt) / T_nd)
        n_cycles_ref = int(ref[-1, 1] / T_nd)

        tr_peaks = per_cycle_peak(tr[:, 1], tr[:, 5], t_ckpt, T_nd, n_cycles_tr)
        ref_peaks = per_cycle_peak(ref[:, 1], ref[:, 5], 0.0, T_nd, n_cycles_ref)

        settle_cycle, reldiff = find_settle(tr_peaks, ref_peaks)
        d_rpm = to["rpm"] - frm["rpm"]
        d_theta = to["theta"] - frm["theta"]
        results.append({
            "level": level, "from": frm, "to": to,
            "d_rpm": d_rpm, "d_theta": d_theta,
            "t_settle": settle_cycle, "n_cycles_available": n_cycles_tr,
        })
        settle_str = str(settle_cycle) if settle_cycle is not None else "NEVER"
        frm_str = f"{frm['rpm']:g}/{frm['theta']:g}"
        to_str = f"{to['rpm']:g}/{to['theta']:g}"
        print(f"{level:>2} {frm_str:>14} {to_str:>14} {d_rpm:>7.1f} {d_theta:>8.1f} {settle_str:>14}")

out_path = ROOT / "experiments" / "settling_study_v2_results.json"
out_path.write_text(json.dumps(results, indent=2))
print(f"\nSaved {out_path}")
