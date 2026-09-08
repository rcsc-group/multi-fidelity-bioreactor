"""The actual Fig 13a chained-sweep readiness test: a real ascending chain,
not independent single hops.

17.5 -> 20 -> 22.5 -> 25 -> 27.5 -> 30 -> 32.5 -> 35 -> 37.5 rpm, theta=7,
L9, video-enabled, each hop warm-started from the PREVIOUS hop's own
result (su per hop 0.867-0.933, matching the natural Fig13a point spacing
-- what a production chained sweep actually does, unlike the independent
large-single-hop test this replaces after user pushback, diary.md
2026-09-08 (5)).

chain.py's build_chain() only supports a single, uniform n_transition_cycles
across every restart segment in the chain (checked directly, not assumed --
per-segment cycle counts aren't part of its interface). N_TRANSITION=25 is
a compromise: generous for the 3 hops being measured against an independent
cold-start reference (22.5, 30, 37.5rpm -- already running/queued as
l9_cold_rpm{22.5,30,37.5}), since these su changes (7-13%) are much milder
than the one L8 data point that needed ~30-35 cycles at su=0.538 (46%); a
bit more than strictly needed for the 5 pass-through hops, which is the
accepted cost of the API not supporting a mixed budget.

8 restart segments * 25 cycles = 200 cycle-equivalents, sequential (each
hop depends on the last) -- at L9's ~11 min/cycle (32 tasks, from
l9_sweep_rpm17.5's known cost) that's roughly 37 hours wall-clock. Uses
chain.py's self-submitting mechanism so it runs unattended; --exclude is
now propagated to every self-submitted segment, not just segment 0 (fixed
alongside this, config/slurm_mpi_template.sh).

Usage:
    uv run python scripts/submit_l9_chain.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain

chain.validate_params = lambda params: None

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-video-fixed"
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
FILL_LEVEL = 0.5
PROJECT_ROOT = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor"
N_TRANSITION = 25

RPMS = [20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5]


def omega_b_of(rpm):
    return rpm * 2 * 3.141592653589793 / 60.0


src_dir = Path(PROJECT_ROOT) / "runs" / "l9_sweep_rpm17.5"
t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])

cfg = {
    "motion": {"omega_b": omega_b_of(RPMS[-1]), "theta_max": [7.0, 0.0, 0.0]},
    "fidelity": 9, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
    "n_mix_cycles": N_TRANSITION,   # unused when initial_checkpoint is set
    # (every segment restarts), but build_chain() requires the key present.
    "n_transition_cycles": N_TRANSITION, "t_buffer": 0.0,
    "sweep": {"parameter": "omega_b", "values": [omega_b_of(r) for r in RPMS]},
    "initial_checkpoint": {
        "t_dump": t_dump, "omega_b": omega_b_of(17.5),
        "theta_max": [7.0, 0.0, 0.0],
        "checkpoint_path": str(src_dir / "checkpoint.dump"),
    },
    "mpi": True, "ntasks": 32, "mem_per_cpu": "4G", "walltime": "10:00:00",
    "binary": BINARY, "submit": True, "exclude": "node2336",
}
run_job_ids = chain.submit_chain(cfg)
print(f"\nL9 chain (rpm order {RPMS}): {run_job_ids}")
print("(each entry is (run_id, job_id); rpm k corresponds to entry k)")
