"""Test C, expanded field list, phase 1 (diary.md 2026-09-07): re-run the
source cold start with the expanded restart_diagnostic instrumentation
(now checks pf, g.x, g.y in addition to u.x/u.y/p/f, which all matched
exactly in the first pass). The kick test ruled out the continuous PDE
itself as the source of the escape-to-a-different-branch behavior, so the
remaining candidates are fields the first Test C pass didn't check.

Usage:
    uv run python scripts/test_C_expanded_fields.py
"""
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
TESTC_BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-testc"

src_params = {
    "run_id": "testc2_source_th4.1", "fidelity": 6,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
    "n_harmonics": 1, "theta_max": [4.1, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": 32.5 * 2 * 3.141592653589793 / 60.0, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 15.183225, "n_mix_cycles": 80, "_binary": TESTC_BINARY,
}
job_id = submit_slurm(src_params, project_root=root, walltime="01:00:00",
                       template=root / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=8, mem="4G", exclude="node2336")
print(f"source (re-run, expanded diagnostic): job={job_id}  run_id=testc2_source_th4.1")
