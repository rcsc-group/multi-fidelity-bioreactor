"""Submit fidelity-6 and fidelity-7 health-check runs to SLURM.

Uses config/slurm_health_template.sh which calls the BioReactor-health binary
(compiled with DIAGNOSTICS=1) to produce pressure_diag.dat in addition to
standard outputs.  Build the health binary first with 'make build-health'.
"""
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
    "t_end": 100.0,
}

root     = Path(__file__).parents[1]
template = root / "config" / "slurm_health_template.sh"

for level, walltime in [(6, "00:45:00"), (7, "02:00:00")]:
    params = {**BASE, "run_id": f"health_l{level}", "fidelity": level}
    job_id = submit_slurm(params, project_root=root, walltime=walltime, template=template)
    print(f"fidelity={level}  job_id={job_id}  walltime={walltime}  run_dir=runs/health_l{level}/")
