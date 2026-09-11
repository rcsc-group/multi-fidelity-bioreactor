"""Rerun with video enabled, to visually inspect the L6->L7 cross-level
pilot's ~31% wrong-plateau finding for a bug. Neither the pilot nor any
existing L6/L7 reference run has video frames (none of them used a
video-enabled binary), so all three are fresh, cheap reruns (minutes at
this fidelity) rather than reusing anything:

  1. crosslevel_L6toL7_pilot_vid -- same as the pilot (L6 source ->
     refine to L7 -> 60 cycles), video-enabled this time.
  2. l7_coldstart_vid -- ordinary L7 cold start, SAME condition and cycle
     count, for a matched-frame-count comparison video.
  3. l6_source_vid -- the L6 source's own dynamics, for context on what's
     being interpolated.

mbessa-condo, extending the user's "for this experiment" authorization to
cover generating the comparison video for the SAME finding, not a new one.

Usage:
    uv run python scripts/submit_cross_level_videos.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-crosslevel-video"
omega = 32.5 * 2 * 3.141592653589793 / 60.0
GEOM = {"a": 0.25, "b": 0.03575, "n": 8.0}

# --- 1. cross-level pilot, with video --------------------------------------
src_dir = root / "runs" / "settling_baseline_L6"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
p1 = {
    "run_id": "crosslevel_L6toL7_pilot_vid", "fidelity": 7,
    "geometry": GEOM, "fill_level": 0.5, "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 60 * 0.607329, "n_mix_cycles": 60, "frames_per_period": 5,
    "t_checkpoint": t_dump,
    "omega_b_prev": omega, "theta_max_prev": [7.0, 0.0, 0.0],
    "phi_angular_prev": [0.0, 0.0, 0.0], "amplitude_h_prev": [0.0, 0.0, 0.0],
    "phi_horizontal_prev": [0.0, 0.0, 0.0], "omega_h_prev": 0.0,
    "_binary": BINARY,
}
j1 = submit_slurm(p1, project_root=root, runs_root=root / "runs",
                  walltime="02:00:00",
                  template=root / "config" / "slurm_mpi_template.sh",
                  checkpoint=str((src_dir / "checkpoint.dump").resolve()),
                  cpus=1, ntasks=8, mem="4G",
                  account="mbessa-condo", qos="mbessa-condo")
print(f"1. crosslevel pilot (video): job={j1}")

# --- 2. matched L7 cold start ------------------------------------------------
p2 = {
    "run_id": "l7_coldstart_vid", "fidelity": 7,
    "geometry": GEOM, "fill_level": 0.5, "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 60 * 0.607329, "n_mix_cycles": 60, "frames_per_period": 5,
    "_binary": BINARY,
}
j2 = submit_slurm(p2, project_root=root, walltime="02:00:00",
                  template=root / "config" / "slurm_mpi_template.sh",
                  cpus=1, ntasks=8, mem="4G",
                  account="mbessa-condo", qos="mbessa-condo")
print(f"2. L7 cold start (video): job={j2}")

# --- 3. L6 source's own video ------------------------------------------------
p3 = {
    "run_id": "l6_source_vid", "fidelity": 6,
    "geometry": GEOM, "fill_level": 0.5, "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega, "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": 60 * 0.607329, "n_mix_cycles": 60, "frames_per_period": 5,
    "_binary": BINARY,
}
j3 = submit_slurm(p3, project_root=root, walltime="02:00:00",
                  template=root / "config" / "slurm_mpi_template.sh",
                  cpus=1, ntasks=8, mem="4G",
                  account="mbessa-condo", qos="mbessa-condo")
print(f"3. L6 source (video): job={j3}")
