"""Resubmit the 6 decreasing-direction dtheta transitions cancelled by
SLURM (reason QOSMaxMemoryPerUser) due to the chain.py mem_per_cpu bug
fixed in this same session (diary.md 2026-09-07) -- theta_source in
{3, 4, 5, 6.5, 6.9, 7}. Updates the existing manifest in place.

Usage:
    uv run python scripts/resubmit_dtheta_decreasing_cancelled.py
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

CANCELLED_THETA_SOURCES = [3.0, 4.0, 5.0, 6.5, 6.9, 7.0]


def source_run_id(theta: float) -> str:
    if theta == 4.0:
        return "settling_ref_L6_rpm32.5_th4"
    if theta == 7.0:
        return "settling_baseline_L6"
    return f"dtheta_src_L6_th{theta:g}"


manifest_path = Path(PROJECT_ROOT) / "experiments" / "dtheta_monotonicity_decreasing_manifest.json"
manifest = json.loads(manifest_path.read_text())

for theta_source in CANCELLED_THETA_SOURCES:
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

manifest_path.write_text(json.dumps(manifest, indent=2))
print(f"\nUpdated {manifest_path}")
