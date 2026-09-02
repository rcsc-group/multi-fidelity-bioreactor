"""Redo of the fig13a_redo 9-point sweep with upstream's ramp matched
(diary.md 2026-09-01, user-flagged regression -- fig13a_redo used our
fork's own 3-cycle smooth-step ramp instead of upstream's 30-second
linear ramp). Same 9 independent cold starts, same L8/t_end=20, only the
ramp mechanism changes (BioReactor_rampmatch_mpi).
"""
import math
import sys

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

PROJECT_ROOT = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor"
TEMPLATE = f"{PROJECT_ROOT}/config/slurm_mpi_template.sh"
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-rampmatch"

RPMS = [17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5]

job_ids = []
for rpm in RPMS:
    omega_b = rpm * 2 * math.pi / 60.0
    run_id = f"fig13a_rampmatch_rpm{rpm:g}"
    params = {
        "run_id": run_id,
        "fidelity": 8,
        "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0},
        "fill_level": 0.5,
        "n_harmonics": 1,
        "theta_max": [7.0, 0.0, 0.0],
        "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": omega_b,
        "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0],
        "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": 20.0,
        "n_mix_cycles": 80,
        "_binary": BINARY,
    }
    job_id = submit_slurm(
        params,
        project_root=PROJECT_ROOT,
        walltime="04:00:00",
        template=TEMPLATE,
        mem="4G",
        cpus=1,
        ntasks=16,
    )
    job_ids.append((run_id, job_id))
    print(f"{run_id}: omega_b={omega_b:.4f} -> job {job_id}")

print()
print("job_ids:", job_ids)
