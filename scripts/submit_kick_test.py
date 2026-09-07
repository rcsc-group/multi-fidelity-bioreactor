"""Velocity-kick test (diary.md 2026-09-07): single, continuous,
uninterrupted cold start at (32.5rpm, theta=7deg), L7 -- NO checkpoint
restart involved at all. The kick-test binary applies a one-shot
u *= (1+1e-4) perturbation at t=KICK_T=12.15 (cycle 20, well after the
flow has settled onto the reference branch), then the run continues
normally for ~120 more cycles. If the trajectory escapes to the same
~1.25x branch the checkpoint-restart cases converge to, that's a genuine
property of the continuous dynamical system, not anything specific to
dump/restore. If it stays on the reference branch, that points back to
something restart-specific.

Usage:
    uv run python scripts/submit_kick_test.py
"""
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, TEMPLATE, omega_b_of

KICK_BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-kicktest"
T_PER_ND = 0.607329
N_CYCLES = 140  # 20 to settle + kick + 120 to observe

params = {
    "run_id": "kicktest_L7_th7",
    "fidelity": 7,
    "geometry": GEOMETRY,
    "fill_level": FILL_LEVEL,
    "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0],
    "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": omega_b_of(32.5),
    "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": N_CYCLES * T_PER_ND,
    "n_mix_cycles": N_CYCLES,
    "_binary": KICK_BINARY,
}
job_id = submit_slurm(params, project_root=Path(PROJECT_ROOT), walltime="02:00:00",
                       template=Path(TEMPLATE), ntasks=8, mem="4G")
print(f"kick test -> run={params['run_id']}  job={job_id}")
