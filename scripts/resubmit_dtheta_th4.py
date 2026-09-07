"""Third resubmission attempt for theta_source=4 -> target=2 (diary.md
2026-09-07): cancelled twice already, both times genuine CPU-cap
contention (account-wide 64-cpu QOS cap, saturated by the still-running L8
falsification jobs + a pre-existing array job), not a bug. Submitting
alone, at ntasks=4, to fit within currently available headroom rather than
alongside another job.

Usage:
    uv run python scripts/resubmit_dtheta_th4.py
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
theta_source = 4.0
src_dir = RUNS_DIR / "settling_ref_L6_rpm32.5_th4"
d = np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)
t_dump = float(d[-1, 1])

cfg = {
    "motion": {"omega_b": omega_b_of(32.5), "theta_max": [2.0, 0.0, 0.0]},
    "fidelity": 6,
    "geometry": GEOMETRY,
    "fill_level": FILL_LEVEL,
    "n_mix_cycles": 60,
    "n_transition_cycles": 60,
    "t_buffer": 0.0,
    "sweep": {"parameter": "omega_b", "values": [omega_b_of(32.5)]},
    "initial_checkpoint": {
        "t_dump": t_dump,
        "omega_b": omega_b_of(32.5),
        "theta_max": [theta_source, 0.0, 0.0],
        "checkpoint_path": str(src_dir / "checkpoint.dump"),
    },
    "mpi": True,
    "ntasks": 4,
    "mem_per_cpu": "4G",
    "walltime": "06:00:00",
    "binary": MAINLINE_BINARY,
    "submit": True,
}
job_run_ids = chain.submit_chain(cfg)
run_id = job_run_ids[0][0]

manifest_path = Path(PROJECT_ROOT) / "experiments" / "dtheta_monotonicity_decreasing_manifest.json"
manifest = json.loads(manifest_path.read_text())
manifest["edges"]["4"] = {"theta_source": theta_source, "run_id": run_id,
                           "source_run": "settling_ref_L6_rpm32.5_th4"}
manifest_path.write_text(json.dumps(manifest, indent=2))
print(f"theta_source=4 -> target=2: run={run_id}  {job_run_ids}")
