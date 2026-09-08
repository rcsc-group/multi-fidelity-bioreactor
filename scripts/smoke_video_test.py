"""Smoke test the video-enabled binary (correctness fixes + VIDEOS=1) at L6
before committing to expensive L9 runs, and measure video overhead.

Usage:
    uv run python scripts/smoke_video_test.py
"""
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
omega = 32.5 * 2 * 3.141592653589793 / 60.0
params = {
    "run_id": "smoke_video_test", "fidelity": 6,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
    "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 5 * 0.607329, "n_mix_cycles": 5, "frames_per_period": 5,
    "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-video-fixed",
}
job = submit_slurm(params, project_root=root, walltime="00:20:00",
                   template=root / "config" / "slurm_mpi_template.sh",
                   cpus=1, ntasks=8, mem="4G", exclude="node2336")
print(f"smoke_video: job={job}")
