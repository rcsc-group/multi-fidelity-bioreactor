"""Does checkpointing (warm start) actually converge faster than a cold
start, on the theta-only edges where there's no frequency-discontinuity
confound? Zero new compute: measures each reference run's OWN cold-start
settling time (its early cycles compared against its own converged
late-cycle value, same rigorous method as the warm-start analysis) and
compares directly against the already-computed warm-start settling times
for the same target condition -- diary.md 2026-09-07.

Usage:
    uv run python scripts/compare_coldstart_vs_warmstart_settling.py
"""
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"
TOL = 0.05

THETA_ONLY_EDGES = [
    ({"rpm": 32.5, "theta": 7.0}, {"rpm": 32.5, "theta": 4.0}),
    ({"rpm": 32.5, "theta": 4.0}, {"rpm": 32.5, "theta": 2.0}),
    ({"rpm": 17.5, "theta": 7.0}, {"rpm": 17.5, "theta": 4.0}),
    ({"rpm": 17.5, "theta": 4.0}, {"rpm": 17.5, "theta": 2.0}),
    ({"rpm": 37.5, "theta": 7.0}, {"rpm": 37.5, "theta": 4.0}),
    ({"rpm": 37.5, "theta": 4.0}, {"rpm": 37.5, "theta": 2.0}),
]


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


def cold_start_settle(rpm, theta, level, tol=TOL, n_tail=5):
    run_id = f"settling_baseline_L{level}" if (rpm, theta) == (32.5, 7.0) else f"settling_ref_L{level}_rpm{rpm:g}_th{theta:g}"
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


v2_manifest = json.loads((ROOT / "experiments" / "settling_study_v2_manifest.json").read_text())
v2_results = {}
for level_str, edges in v2_manifest["edges"].items():
    for tag, e in edges.items():
        key = (level_str, e["from"]["rpm"], e["from"]["theta"], e["to"]["rpm"], e["to"]["theta"])
        v2_results[key] = e["run_id"]

# rerun v2 settle values directly (avoid depending on a separate results file being current)
import importlib.util
spec = importlib.util.spec_from_file_location("analyze_v2", ROOT / "scripts" / "analyze_settling_study_v2.py")
# instead of importing (which executes prints), just recompute inline using the same logic
def warmstart_settle(frm, to, level, run_id, tol=TOL, n_tail=5):
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
    target = float(np.nanmean(ref_peaks[-n_tail:]))
    reldiff = np.abs(tr_peaks - target) / target
    for start in range(len(reldiff)):
        if np.all(reldiff[start:] < tol):
            return start
    return None


print(f"{'L':>2} {'from':>10} {'to':>10} {'coldstart':>10} {'warmstart':>10} {'faster?':>8}")
for level in [6, 7, 8]:
    for frm, to in THETA_ONLY_EDGES:
        run_id = v2_manifest["edges"][str(level)][
            f"{frm['rpm']:g}_{frm['theta']:g}__{to['rpm']:g}_{to['theta']:g}"
        ]["run_id"]
        cold = cold_start_settle(to["rpm"], to["theta"], level)
        warm = warmstart_settle(frm, to, level, run_id)
        cold_str = str(cold) if cold is not None else "NEVER"
        warm_str = str(warm) if warm is not None else "NEVER"
        if cold is not None and warm is not None:
            faster = "YES" if warm < cold else ("SAME" if warm == cold else "NO")
        else:
            faster = "?"
        frm_str = f"{frm['rpm']:g}/{frm['theta']:g}"
        to_str = f"{to['rpm']:g}/{to['theta']:g}"
        print(f"{level:>2} {frm_str:>10} {to_str:>10} {cold_str:>10} {warm_str:>10} {faster:>8}")
