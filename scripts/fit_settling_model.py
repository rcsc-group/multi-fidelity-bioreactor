"""Build the consolidated settling model fit: adaptive-tolerance settle
cycle for every measured transition, using extended (up to 180-cycle) data
where we ran it (scripts/extend_censored_settling_runs*.py) and the
original 60-cycle data everywhere else.

Excludes 32.5rpm theta 7->2 at ALL levels: this specific transition,
extended to 180 cycles, plateaued at a level distinctly above its own
reference's noise floor (~2x) instead of continuing to converge -- flagged
2026-09-12 as a likely genuine distinct-branch anomaly (same family as the
cross-level warm-start finding earlier this session), not a settling-time
question. Do not feed it into a cycle-count model; it needs its own
investigation.

Usage:
    uv run python scripts/fit_settling_model.py
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1]))
from scripts.settling_model import (
    t_per_nd, waveforms_from_series, residual_series, adaptive_tol, _settle_from_residual,
)

ROOT = Path(__file__).parent.parent
RUNS = ROOT / "runs"

# Conditions extended past the original 60 cycles -- (level, rpm, theta) -> [run_id chain]
EXTENDED_CHAINS = {
    (6, 32.5, 4.0): ["7bb073d2", "fd00e7b8", "b5e64ccb"],   # 180 cyc, converged
    (8, 32.5, 4.0): ["865f5441", "c35b6eee"],                # 120 cyc, converged (~91)
    (7, 32.5, 2.0): ["2c107144", "e2a0fe1e", "c420a2b5"],   # 180 cyc, FLAGGED -- excluded below
    (6, 17.5, 4.0): ["b31a1248", "f081ff04"],                # 120 cyc, converged (matches ref floor)
    (8, 37.5, 7.0): ["77a88c02", "354808af"],                # 120 cyc, converged (matches ref floor)
    (8, 37.5, 2.0): ["70db6266", "f0839974"],                # 120 cyc, borderline-converged
    (7, 17.5, 2.0): ["eb021b15", "1a899cc2"],                # 120 cyc, converged (matches ref floor)
}

FLAGGED_ANOMALY = {(32.5, 4.0, 32.5, 2.0)}  # (from_rpm, from_theta, to_rpm, to_theta) -- all levels
# NOTE: the anomalous edge is 32.5rpm theta 4->2 (one hop of the chained
# baseline(7) -> 4 -> 2 topology), NOT a direct 7->2 transition -- the
# manifest never has that edge. Caught a mislabeled first version of this
# key (7.0 instead of 4.0) by testing settling_cycles() against it directly
# and finding it silently failed to raise.


def _load_series(run_ids):
    ts, ys = [], []
    for rid in run_ids:
        d = np.loadtxt(RUNS / rid / "shear_stress.dat", skiprows=1)
        ts.append(d[:, 1])
        ys.append(d[:, 5])
    return np.concatenate(ts), np.concatenate(ys)


def fit() -> list[dict]:
    manifest = json.loads((ROOT / "experiments" / "settling_study_v2_manifest.json").read_text())
    results = []
    for level_str, edges in manifest["edges"].items():
        level = int(level_str)
        for tag, e in edges.items():
            frm, to, run_id = e["from"], e["to"], e["run_id"]
            key = (frm["rpm"], frm["theta"], to["rpm"], to["theta"])
            if key in FLAGGED_ANOMALY:
                results.append({
                    "level": level, "from": frm, "to": to,
                    "d_rpm": to["rpm"] - frm["rpm"], "d_theta": to["theta"] - frm["theta"],
                    "settle_cycle": None, "confidence": "flagged_anomaly",
                    "note": "excluded: plateaus above reference noise floor after 180 cycles, "
                            "see diary 2026-09-12 -- do not use for cycle-count sizing",
                })
                continue

            ext_key = (level, to["rpm"], to["theta"])
            if ext_key in EXTENDED_CHAINS:
                run_ids = EXTENDED_CHAINS[ext_key]
                t, y = _load_series(run_ids)
            else:
                d = np.loadtxt(RUNS / run_id / "shear_stress.dat", skiprows=1)
                t, y = d[:, 1], d[:, 5]

            params = json.loads((RUNS / run_id / "params.json").read_text())
            t_ckpt = params["t_checkpoint"]
            ref_id = f"settling_ref_L{level}_rpm{to['rpm']:g}_th{to['theta']:g}"
            ref_dir = RUNS / ref_id
            tol = adaptive_tol(ref_dir, to["rpm"], to["theta"])
            resid, _ = residual_series(t, y, ref_dir, to["rpm"], to["theta"], t_checkpoint=t_ckpt)
            k = _settle_from_residual(resid, tol)
            n_avail = len(resid) if resid is not None else 0

            results.append({
                "level": level, "from": frm, "to": to,
                "d_rpm": to["rpm"] - frm["rpm"], "d_theta": to["theta"] - frm["theta"],
                "settle_cycle": k, "n_cycles_available": n_avail, "tol_used": tol,
                "confidence": "measured" if k is not None else "censored_never_settled",
            })
    return results


if __name__ == "__main__":
    rows = fit()
    print(f"{'L':>2} {'from':>14} {'to':>14} {'settle_cycle':>12} {'tol':>6} {'confidence':>20}")
    for r in rows:
        frm_s = f"{r['from']['rpm']:g}/{r['from']['theta']:g}"
        to_s = f"{r['to']['rpm']:g}/{r['to']['theta']:g}"
        sc = "NEVER" if r["settle_cycle"] is None else r["settle_cycle"]
        tol = f"{r.get('tol_used', float('nan')):.1%}" if r.get("tol_used") else "--"
        print(f"{r['level']:>2} {frm_s:>14} {to_s:>14} {sc!s:>12} {tol:>6} {r['confidence']:>20}")

    out = Path(__file__).parent.parent / "experiments" / "settling_model_fit.json"
    out.write_text(json.dumps(rows, indent=2))
    print(f"\nsaved {out}")

    measured = [r for r in rows if r["confidence"] == "measured"]
    censored = [r for r in rows if r["confidence"] == "censored_never_settled"]
    print(f"\n{len(measured)} measured, {len(censored)} still censored (never settled even with "
          f"extension, excluding the flagged anomaly), max measured settle_cycle="
          f"{max((r['settle_cycle'] for r in measured), default=None)}")
