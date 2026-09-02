"""L6 arm of the Fig 13a replica (diary.md 2026-09-02): same 9-point,
theta=7deg sweep as fig13a_rampmatch, same ramp-matched current driver
(BioReactor-mpi-rampmatch, fidelity is a runtime param -- 1<<fidelity
grid cells/side), just fidelity=6 instead of 8. 9 independent cold
starts, not sweep.py chaining, for the same restart-transient-avoidance
reason as the L8 run.

Submitted on mbessa-condo per user's explicit per-job permission
(2026-09-02). ntasks=4 (not 16) -- L6 is a small enough grid that 16
MPI ranks would mostly add domain-decomposition overhead, not speed.
"""
import math
import subprocess
import sys

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

PROJECT_ROOT = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor"
TEMPLATE = f"{PROJECT_ROOT}/config/slurm_mpi_template.sh"
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-rampmatch"
SCRATCH_RUNS = "/oscar/scratch/eaguerov/mpi_runs"

import os
RPMS = [32.5] if os.environ.get("SMOKE_TEST") else [17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5]

job_ids = []
for rpm in RPMS:
    omega_b = rpm * 2 * math.pi / 60.0
    run_id = f"fig13a_l6_rpm{rpm:g}"
    params = {
        "run_id": run_id,
        "fidelity": 6,
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
    # Stage via submit_slurm() (normal QOS) so scratch params.json / canonical
    # run dir get set up correctly, then cancel and resubmit on mbessa-condo --
    # same two-step pattern as fig13a_rampmatch (submit_slurm has no
    # account/QOS override). Re-pass the same resource flags submit_slurm()
    # used (it applies them as sbatch CLI flags, not baked into the template
    # -- a bare account/qos-only resubmit would silently fall back to the
    # template's own #SBATCH defaults, ntasks=16/time=04:00:00).
    walltime, mem, cpus, ntasks = "01:00:00", "2G", 1, 4
    job_id = submit_slurm(
        params,
        project_root=PROJECT_ROOT,
        walltime=walltime,
        template=TEMPLATE,
        mem=mem,
        cpus=cpus,
        ntasks=ntasks,
    )
    subprocess.run(["scancel", job_id], check=True)

    params_path = f"{SCRATCH_RUNS}/{run_id}/params.json"
    result = subprocess.run(
        [
            "sbatch", "--account=mbessa-condo", "--qos=mbessa-condo",
            "--parsable", "--no-requeue",
            f"--time={walltime}",
            f"--mem-per-cpu={mem}",
            f"--cpus-per-task={cpus}",
            f"--ntasks={ntasks}",
            f"--export=NONE,PARAMS={params_path}",
            TEMPLATE,
        ],
        cwd=PROJECT_ROOT, capture_output=True, text=True, check=True,
    )
    condo_job_id = result.stdout.strip()
    job_ids.append((run_id, condo_job_id))
    print(f"{run_id}: omega_b={omega_b:.4f} -> job {condo_job_id} (mbessa-condo)")

print()
print("job_ids:", job_ids)
