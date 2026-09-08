"""Fig 13a readiness check at L9: 3 targets, warm (single-hop restart from
the existing 17.5rpm L9 source) vs cold (independent start), video-enabled.

Targets: 22.5, 30.0, 37.5 rpm (theta=7, L9) -- su from 17.5rpm source =
0.778, 0.583, 0.467, all comfortably above the L7 escape threshold
(0.126-0.154) and squarely in the range the L8 check (su=0.538) showed
converges cleanly. This extends that single L8 data point to 3 more at L9.

Source checkpoint (l9_sweep_rpm17.5) was dumped by a pre-fix, non-video
binary; verified safe to restore into the video-enabled, fully-fixed binary
-- Basilisk's restore() matches fields by name from the file, so the extra
video-only scalars (tau_field, ediss_field) are simply left untouched, not
misaligned.

n_mix_cycles=45 (video-fixed, phase-fixed, g-mode-2 binary), ntasks=32
(matches the known-good L9 config from l9_sweep_rpm17.5), frames_per_period=5.

Usage:
    uv run python scripts/submit_l9_warmcold_compare.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-video-fixed"
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
FILL_LEVEL = 0.5
T_PER_ND = 0.607329
N_MIX = 45


def omega_b_of(rpm):
    return rpm * 2 * 3.141592653589793 / 60.0


TARGETS = [22.5, 30.0, 37.5]
src_dir = root / "runs" / "l9_sweep_rpm17.5"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])

jobs = {}
for rpm in TARGETS:
    # cold start
    cold_params = {
        "run_id": f"l9_cold_rpm{rpm:g}", "fidelity": 9,
        "geometry": GEOMETRY, "fill_level": FILL_LEVEL, "n_harmonics": 1,
        "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": omega_b_of(rpm), "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": N_MIX * T_PER_ND, "n_mix_cycles": N_MIX, "frames_per_period": 5,
        "_binary": BINARY,
    }
    j = submit_slurm(cold_params, project_root=root, walltime="18:00:00",
                     template=root / "config" / "slurm_mpi_template.sh",
                     cpus=1, ntasks=32, mem="4G", exclude="node2336")
    jobs[f"cold_{rpm:g}"] = j
    print(f"cold {rpm:g}rpm: job={j}  run_id=l9_cold_rpm{rpm:g}")

    # warm restart, single hop from 17.5rpm source
    warm_params = {
        "run_id": f"l9_warm_rpm{rpm:g}", "fidelity": 9,
        "geometry": GEOMETRY, "fill_level": FILL_LEVEL, "n_harmonics": 1,
        "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": omega_b_of(rpm), "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": N_MIX * T_PER_ND, "n_mix_cycles": N_MIX, "frames_per_period": 5,
        "t_checkpoint": t_dump,
        "omega_b_prev": omega_b_of(17.5), "theta_max_prev": [7.0, 0.0, 0.0],
        "phi_angular_prev": [0.0, 0.0, 0.0], "amplitude_h_prev": [0.0, 0.0, 0.0],
        "phi_horizontal_prev": [0.0, 0.0, 0.0], "omega_h_prev": 0.0,
        "_binary": BINARY,
    }
    j = submit_slurm(warm_params, project_root=root, runs_root=root / "runs",
                     walltime="18:00:00",
                     template=root / "config" / "slurm_mpi_template.sh",
                     checkpoint=str((src_dir / "checkpoint.dump").resolve()),
                     cpus=1, ntasks=32, mem="4G", exclude="node2336")
    jobs[f"warm_{rpm:g}"] = j
    print(f"warm {rpm:g}rpm (su={17.5/rpm:.3f}): job={j}  run_id=l9_warm_rpm{rpm:g}")

print("\nall jobs:", jobs)
