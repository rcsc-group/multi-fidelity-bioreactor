"""Recompute settle_cycle for the 7 representative censored transitions
after extending each by another 60 cycles from its own checkpoint (see
scripts/extend_censored_settling_runs.py) -- concatenates the source run's
shear_stress.dat with its extension's (same absolute time axis, since the
extension restarts exactly where the source left off) and re-runs the
validated waveform-RMS settle detection over the full ~120-cycle window.

Usage:
    uv run python scripts/analyze_extended_settling.py
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1]))
from scripts.settling_model import (
    t_per_nd, waveforms_from_series, converged_waveform, best_phase_shift, NPH,
)

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"


def settle_cycle_extended(source_id, ext_id, level, rpm_to, theta_to, t_checkpoint,
                           tol, n_tail=20, nph=NPH):
    d_src = np.loadtxt(RUNS / source_id / "shear_stress.dat", skiprows=1)
    d_ext = np.loadtxt(RUNS / ext_id / "shear_stress.dat", skiprows=1)
    t = np.concatenate([d_src[:, 1], d_ext[:, 1]])
    y = np.concatenate([d_src[:, 5], d_ext[:, 5]])

    T = t_per_nd(rpm_to, theta_to)
    ref_id = f"settling_ref_L{level}_rpm{rpm_to:g}_th{theta_to:g}"
    ref_waves = waveforms_from_series(*np.loadtxt(RUNS / ref_id / "shear_stress.dat", skiprows=1)[:, [1, 5]].T, T)
    ref_conv = converged_waveform(ref_waves, n_tail)
    scale = np.sqrt(np.mean(ref_conv ** 2))

    run_waves = waveforms_from_series(t, y, T, t0=t_checkpoint, nph=nph)
    run_conv = converged_waveform(run_waves, n_tail)
    shift = best_phase_shift(run_conv, ref_conv, nph)
    phase = np.linspace(0, 1, nph, endpoint=False)
    q = np.mod(phase + shift, 1.0)
    resid = np.full(len(run_waves), np.nan)
    for k, w in enumerate(run_waves):
        if w is None:
            continue
        w_aligned = np.interp(q, phase, w, period=1.0)
        resid[k] = np.sqrt(np.mean((w_aligned - ref_conv) ** 2)) / scale

    for k in range(len(resid)):
        seg = resid[k:]
        valid = seg[~np.isnan(seg)]
        if len(valid) > 0 and np.all(valid < tol):
            return k, len(run_waves)
    return None, len(run_waves)


manifest = json.loads((ROOT / "experiments" / "settling_extension_manifest.json").read_text())

print(f"{'L':>2} {'condition':>14} {'tol':>5} {'settle_cycle':>12} {'n_avail':>8}")
results = []
for e in manifest:
    level = e["level"]
    source_id, ext_id = e["source_run_id"], e["extension_run_id"]
    p_src = json.loads((RUNS / source_id / "params.json").read_text())
    t_ck = p_src["t_checkpoint"]
    rpm, theta = e["rpm"], e["theta"]
    for tol in (0.05, 0.02, 0.01):
        k, n = settle_cycle_extended(source_id, ext_id, level, rpm, theta, t_ck, tol)
        sc = "NEVER" if k is None else k
        print(f"{level:>2} {f'{rpm:g}/{theta:g}':>14} {tol:>5.0%} {sc!s:>12} {n:>8}")
        results.append({"level": level, "rpm": rpm, "theta": theta, "tol": tol,
                         "settle_cycle": k, "n_cycles_available": n,
                         "source_run_id": source_id, "extension_run_id": ext_id})

out = ROOT / "experiments" / "settling_extended_results.json"
out.write_text(json.dumps(results, indent=2))
print(f"\nsaved {out}")
