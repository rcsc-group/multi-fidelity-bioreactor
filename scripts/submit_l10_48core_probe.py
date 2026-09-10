"""48-core (single node) L10 throughput calibration before committing to a
9-point Fig13a sweep at this exact configuration. The 240-task/5-node probe
already showed multi-node output corruption (diary.md 2026-09-10) -- 48
cores stays on ONE node (this cluster's nodes have 48 cores each),
avoiding that risk entirely, but has never been measured at L10.

5 cycles, cold start, 32.5rpm/theta=7 (same condition as the 240-task
probe, for a clean comparison).

Usage:
    uv run python scripts/submit_l10_48core_probe.py
"""
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
omega = 32.5 * 2 * 3.141592653589793 / 60.0
params = {
    "run_id": "l10_48core_probe_rpm32.5", "fidelity": 10,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
    "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 5 * 0.607329, "n_mix_cycles": 5, "frames_per_period": 5,
    "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-video-fixed",
}
job = submit_slurm(params, project_root=root, walltime="06:00:00",
                   template=root / "config" / "slurm_mpi_template.sh",
                   cpus=1, ntasks=48, mem="6G", exclude="node2336")
print(f"48-core L10 probe: job={job}")
