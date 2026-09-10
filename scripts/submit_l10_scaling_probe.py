"""L10 scaling probe + real progress, combined.

Never run this problem past 64 cores -- every historical L10 job used 32
or 64. "Finish an L10 point in a day at 256+ cores" would be pure
extrapolation. This measures real throughput at high core count directly
from Basilisk's own end-of-run print (points.step/s), while also making
real progress on an actually-useful point: 32.5rpm/theta=7, Kim's
canonical baseline, matching existing data at every other fidelity level.

10 cycles, cold start, 256 tasks (leaves ~56 of the new 312-CPU cap free
for the theta-sweep chain to run alongside it). Video-enabled, fully-fixed
binary. Will be checked by scripts/validate_run.py before any number from
it is reported, same as everything else.

Usage:
    uv run python scripts/submit_l10_scaling_probe.py
"""
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
omega = 32.5 * 2 * 3.141592653589793 / 60.0
params = {
    "run_id": "l10_scaling_probe_rpm32.5", "fidelity": 10,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
    "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 10 * 0.607329, "n_mix_cycles": 10, "frames_per_period": 5,
    "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-video-fixed",
}
job = submit_slurm(params, project_root=root, walltime="06:00:00",
                   template=root / "config" / "slurm_mpi_template.sh",
                   cpus=1, ntasks=256, mem="6G", exclude="node2336")
print(f"L10 scaling probe: job={job}")
