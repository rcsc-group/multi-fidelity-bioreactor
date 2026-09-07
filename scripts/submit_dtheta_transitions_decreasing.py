"""Delta-theta_max monotonicity study, DECREASING-direction mirror, phase 2
(diary.md 2026-09-07): warm-start transitions from theta_source > 2 down
into a fixed target (32.5rpm, theta=2deg, L6) -- the mirror of
submit_dtheta_transitions.py's increasing-direction sweep into theta=7.
Tests whether the settling-time-vs-Delta_theta relationship is direction-
dependent (theta-decreasing is the direction the ORIGINAL catastrophic
anomaly was found in: 7->4->2).

Run submit_dtheta_sources_decreasing.py first and wait for it to complete
(only 2.1, 2.5 are new; 3,4,5,6,6.5,6.9,7 reuse existing converged L6
references from the increasing-direction sweep / parent grid study).

theta_source=2.0 is the trivial Delta_theta_max=0 case (target's own
checkpoint into itself).

Usage:
    uv run python scripts/submit_dtheta_transitions_decreasing.py
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import PROJECT_ROOT, MAINLINE_BINARY, GEOMETRY, FILL_LEVEL, omega_b_of

chain.validate_params = lambda params: None

RUNS_DIR = Path(PROJECT_ROOT) / "runs"
RPM = 32.5
TARGET_THETA = 2.0
FIDELITY = 6
NTASKS = 8
CYCLES_PER_RUN = 60

THETA_SOURCES = [2.0, 2.1, 2.5, 3.0, 4.0, 5.0, 6.0, 6.5, 6.9, 7.0]


def source_run_id(theta: float) -> str:
    if theta == 2.0:
        return "settling_ref_L6_rpm32.5_th2"
    if theta == 4.0:
        return "settling_ref_L6_rpm32.5_th4"
    if theta == 7.0:
        return "settling_baseline_L6"
    return f"dtheta_src_L6_th{theta:g}"


manifest = {"target_theta": TARGET_THETA, "rpm": RPM, "fidelity": FIDELITY, "edges": {}}

for theta_source in THETA_SOURCES:
    src_dir = RUNS_DIR / source_run_id(theta_source)
    d = np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)
    t_dump = float(d[-1, 1])

    cfg = {
        "motion": {"omega_b": omega_b_of(RPM), "theta_max": [TARGET_THETA, 0.0, 0.0]},
        "fidelity": FIDELITY,
        "geometry": GEOMETRY,
        "fill_level": FILL_LEVEL,
        "n_mix_cycles": CYCLES_PER_RUN,
        "n_transition_cycles": CYCLES_PER_RUN,
        "t_buffer": 0.0,
        "sweep": {"parameter": "omega_b", "values": [omega_b_of(RPM)]},
        "initial_checkpoint": {
            "t_dump": t_dump,
            "omega_b": omega_b_of(RPM),
            "theta_max": [theta_source, 0.0, 0.0],
            "checkpoint_path": str(src_dir / "checkpoint.dump"),
        },
        "mpi": True,
        "ntasks": NTASKS,
        "mem_per_cpu": "4G",
        "walltime": "06:00:00",
        "binary": MAINLINE_BINARY,
        "submit": True,
    }
    job_run_ids = chain.submit_chain(cfg)
    run_id = job_run_ids[0][0]
    manifest["edges"][f"{theta_source:g}"] = {
        "theta_source": theta_source, "run_id": run_id, "source_run": source_run_id(theta_source),
    }
    print(f"theta_source={theta_source:g} -> target={TARGET_THETA:g}: run={run_id}  {job_run_ids}")

out_path = Path(PROJECT_ROOT) / "experiments" / "dtheta_monotonicity_decreasing_manifest.json"
out_path.write_text(json.dumps(manifest, indent=2))
print(f"\nSaved {out_path}")
