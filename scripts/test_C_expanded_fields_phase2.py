"""Test C, expanded field list, phase 2 (diary.md 2026-09-07): the two
restarts (trivial and anomalous), from testc2_source_th4.1's checkpoint.
Run test_C_expanded_fields.py first and wait for it to complete.

Usage:
    uv run python scripts/test_C_expanded_fields_phase2.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
src_dir = root / "runs" / "testc2_source_th4.1"
TESTC_BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-testc"

d = np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)
t_dump = float(d[-1, 1])

for tag, target_theta in [("trivial", 4.1), ("anomalous", 4.05)]:
    params = {
        "run_id": f"testc2_{tag}", "fidelity": 6,
        "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
        "n_harmonics": 1, "theta_max": [target_theta, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": 32.5 * 2 * 3.141592653589793 / 60.0, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": 0.6, "n_mix_cycles": 1,
        "t_checkpoint": t_dump,
        "omega_b_prev": 32.5 * 2 * 3.141592653589793 / 60.0,
        "theta_max_prev": [4.1, 0.0, 0.0], "phi_angular_prev": [0.0, 0.0, 0.0],
        "amplitude_h_prev": [0.0, 0.0, 0.0], "phi_horizontal_prev": [0.0, 0.0, 0.0],
        "omega_h_prev": 0.0, "_binary": TESTC_BINARY,
    }
    job_id = submit_slurm(
        params, project_root=root, runs_root=root / "runs", walltime="01:00:00",
        template=root / "config" / "slurm_mpi_template.sh",
        checkpoint=str((src_dir / "checkpoint.dump").resolve()),
        cpus=1, ntasks=8, mem="4G", exclude="node2336",
    )
    print(f"{tag} (target={target_theta:g}): job={job_id}  run_id=testc2_{tag}")
