"""Submit fidelity-6 and fidelity-7 health-check runs to SLURM."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))
from scripts.simulate import submit_slurm

BASE = {
    "omega_b": 3.93, "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_h": 0.0, "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
    "fill_level": 0.5,
    "t_end": 20.0,
}

root = Path(__file__).parents[1]

for level, walltime in [(6, "02:00:00"), (7, "10:00:00")]:
    params = {**BASE, "run_id": f"health_l{level}", "fidelity": level}
    job_id = submit_slurm(params, project_root=root, walltime=walltime)
    print(f"fidelity={level}  job_id={job_id}  walltime={walltime}  run_dir=runs/health_l{level}/")
