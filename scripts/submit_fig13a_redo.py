"""Redo of the Fig 13a(a) replica (tau/EDR vs RPM at theta=7deg, L8),
using the current bug-fixed driver (post tau-histogram race fix, post
Fig 8a sign fix, post H_bio fix) -- diary.md 2026-09-01.

9 INDEPENDENT cold-start runs (not sweep.py's checkpoint-chaining, which
would warm-start each RPM point from the previous one's end state --
exactly the restart-transient confound this session spent real effort
characterizing). t_end=20.0 non-dim (~33 cycles, past Kim et al.'s own
~30-cycle convergence threshold -- an improvement over whatever window
the original Aug 3-4 replica used).
"""
import math
import sys

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

PROJECT_ROOT = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor"
TEMPLATE = f"{PROJECT_ROOT}/config/slurm_mpi_template.sh"
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-fig13a-redo"

RPMS = [17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5]

job_ids = []
for rpm in RPMS:
    omega_b = rpm * 2 * math.pi / 60.0
    run_id = f"fig13a_redo_rpm{rpm:g}"
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
        "n_mix_cycles": 80,  # oxygen/tracer never actually starts (t_end < t_mix) -- this run is for tau_*_qss/max only
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
