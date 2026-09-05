"""L9 sweep, condition 2 of 9 (20rpm), CHECKPOINTED (warm-started) from
condition 1's (17.5rpm) converged end state -- diary.md 2026-09-05.

Uses the MAINLINE binary (BioReactor-mpi), NOT the ramp-matched one
condition 1 used: the ramp-matched acceleration() event hardcodes the
CURRENT omega_b directly (w_bio_st) with no reference to omega_b_prev --
warm-starting with it would jump the forcing frequency instantaneously
at the restart instant. The mainline binary's own N_RAMP_CYCLES=3
smooth-step interpolates both amplitude and phase between the previous
and current omega_b, which is the mechanism this checkpointing step
actually needs. Kim's own methodology never does cross-condition
warm-starting at all (every point is an independent cold start), so
there's no "match Kim" constraint on this choice.

Cycle budget (derived, not guessed -- see diary.md 2026-09-05):
  - N_RAMP_CYCLES=3 (mechanical forcing interpolation, driver default)
  - ~8-9 cycle settling (steady-streaming slow-time argument: 1/theta_max_rad
    = 8.19 cycles for theta_max=7deg -- matches the 2026-08-20 finding's
    empirically-measured ~8-9 cycles for a theta-change warm-start;
    theta_max is unchanged here, so the same mechanism should apply),
    padded to 10 for margin
  - ~25-30 QSS sampling cycles, matching the cold-start convention
  Total: n_transition_cycles=40.

At the measured 21.2 min/cycle (32 ranks, from condition 1): ~14.1h,
comfortable margin under the 48h Exploratory-tier cap.

Usage:
    uv run python scripts/submit_l9_rpm20_checkpoint.py
"""
import math
import sys

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain

# Same validate_params bypass as every other Kim-condition script this
# session -- Kim's literal geometry sits outside config/param_space.yaml's
# BO search bounds by design.
chain.validate_params = lambda params: None

PROJECT_ROOT = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor"
TEMPLATE = f"{PROJECT_ROOT}/config/slurm_mpi_template.sh"
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi"  # mainline, not rampmatched

omega_b_17_5 = 17.5 * 2 * math.pi / 60.0
omega_b_20 = 20.0 * 2 * math.pi / 60.0

cfg = {
    "motion": {"omega_b": omega_b_20, "theta_max": [7.0, 0.0, 0.0]},
    "fidelity": 9,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0},
    "fill_level": 0.5,
    "n_mix_cycles": 40,          # unused when initial_checkpoint is set (n_transition_cycles governs)
    "n_transition_cycles": 40,   # 3 (mechanical) + ~10 (settling, padded) + ~27 (QSS sampling)
    "t_buffer": 0.0,
    "sweep": {"parameter": "omega_b", "values": [omega_b_20]},
    "initial_checkpoint": {
        "t_dump": 24.9,  # exact final t of l9_sweep_rpm17.5 (job 5827446)
        "omega_b": omega_b_17_5,
        "theta_max": [7.0, 0.0, 0.0],
        "checkpoint_path": f"{PROJECT_ROOT}/runs/l9_sweep_rpm17.5/checkpoint.dump",
    },
    "mpi": True,
    "ntasks": 32,
    "mem_per_cpu": "4G",
    "walltime": "24:00:00",
    "binary": BINARY,
    "submit": True,
}

# build_chain() assigns run_id via uuid4 internally -- can't rename it
# through submit_chain()'s cfg dict, so just capture whatever comes back
# and record it clearly here instead of fighting the API.
job_run_ids = chain.submit_chain(cfg)
print(f"L9 condition 2 (20rpm, checkpointed from 17.5rpm): {job_run_ids}")
print(f"  -> canonical run dir: runs/{job_run_ids[0][0]}/")
