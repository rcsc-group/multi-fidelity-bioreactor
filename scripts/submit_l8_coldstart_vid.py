"""Matched L8 cold-start video reference (settling_baseline_L8 predates
video support). Same condition, same 60-cycle window as
crosslevel_L7toL8_pilot, for the two-basin evidence video.

Usage:
    uv run python scripts/submit_l8_coldstart_vid.py
"""
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
omega = 32.5 * 2 * 3.141592653589793 / 60.0
params = {
    "run_id": "l8_coldstart_vid", "fidelity": 8,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
    "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 60 * 0.607329, "n_mix_cycles": 60, "frames_per_period": 5,
    "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-crosslevel-video",
}
job = submit_slurm(params, project_root=root, walltime="02:00:00",
                   template=root / "config" / "slurm_mpi_template.sh",
                   cpus=1, ntasks=8, mem="4G",
                   account="mbessa-condo", qos="mbessa-condo")
print(f"L8 cold start (video): job={job}")
