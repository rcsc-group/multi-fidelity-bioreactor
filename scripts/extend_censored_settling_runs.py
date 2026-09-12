"""Extend a representative subset of the settling-study grid's
right-censored transitions (never settled within their recorded 60 cycles
per scripts/settling_model.py's corrected metric) by another 60 cycles,
continuing from each run's own checkpoint at the SAME condition (a trivial,
already-validated restart -- see "chain seg1 (trivial)" in
scripts/settling_model.py's regression check against diary.md 2026-09-08).

Picked to span levels and cover the most interesting cases:
  - L6/L8 32.5->32.5rpm, theta 7->4  (never settled at ANY of L6/L7/L8)
  - L7    32.5->32.5rpm, theta 4->2  (never settled at ANY level)
  - L6    17.5rpm, theta 7->4        (never at L6, but settled at L7/L8 --
                                       check if L6 just needs longer)
  - L8    32.5->37.5rpm, theta 7     (settles almost instantly at L6/L7 but
                                       NEVER at L8 -- same flavor as this
                                       session's L7->L8 cross-level anomaly)
  - L8    37.5rpm, theta 4->2        (never at ANY level)
  - L7    17.5rpm, theta 4->2        (never at L6/L7, settled at L8)

Walltime sized from each candidate's OWN measured elapsed time for its first
60 cycles (sacct, not extrapolated), with a >=2.5x safety margin -- these are
same-condition continuations, so cost/cycle should closely match the first
segment; margin covers node-contention variance (2-4x observed this session
at other conditions), not model uncertainty.

Usage:
    uv run python scripts/extend_censored_settling_runs.py
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, MAINLINE_BINARY, omega_b_of

chain.validate_params = lambda params: None

PROJECT_ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
RUNS = PROJECT_ROOT / "runs"

# (level, run_id, rpm, theta, walltime) -- walltime from each run's own
# sacct Elapsed for its first 60 cycles, x2.5-3 margin, rounded up.
CANDIDATES = [
    (6, "7bb073d2", 32.5, 4.0, "00:30:00"),   # elapsed 00:02:34
    (8, "865f5441", 32.5, 4.0, "05:00:00"),   # elapsed 01:32:09
    (7, "2c107144", 32.5, 2.0, "01:00:00"),   # elapsed 00:10:11
    (6, "b31a1248", 17.5, 4.0, "00:30:00"),   # elapsed 00:04:14
    (8, "77a88c02", 37.5, 7.0, "05:00:00"),   # elapsed 01:31:37
    (8, "70db6266", 37.5, 2.0, "05:00:00"),   # elapsed 01:53:53
    (7, "eb021b15", 17.5, 2.0, "01:00:00"),   # elapsed 00:13:38
]

EXTEND_CYCLES = 60

results = []
for level, run_id, rpm, theta, walltime in CANDIDATES:
    run_dir = RUNS / run_id
    p = json.loads((run_dir / "params.json").read_text())
    ckpt = run_dir / "checkpoint.dump"
    if not ckpt.exists():
        print(f"SKIP L{level} {run_id}: no checkpoint.dump")
        continue

    t_dump = p["t_checkpoint"] + p["t_end"]  # absolute time this checkpoint was written at

    cfg = {
        "motion": {"omega_b": omega_b_of(rpm), "theta_max": [theta, 0.0, 0.0]},
        "fidelity": level,
        "geometry": GEOMETRY,
        "fill_level": FILL_LEVEL,
        "n_mix_cycles": EXTEND_CYCLES,       # unused when initial_checkpoint is set
        "n_transition_cycles": EXTEND_CYCLES,
        "t_buffer": 0.0,
        "sweep": {"parameter": "omega_b", "values": [omega_b_of(rpm)]},  # trivial: no change
        "initial_checkpoint": {
            "t_dump": t_dump,
            "omega_b": omega_b_of(rpm),
            "theta_max": [theta, 0.0, 0.0],
            "checkpoint_path": str(ckpt.resolve()),
        },
        "mpi": True,
        "ntasks": 8,
        "mem_per_cpu": "4G",
        "walltime": walltime,
        "binary": MAINLINE_BINARY,
        "submit": True,
    }
    job_run_ids = chain.submit_chain(cfg)
    new_run_id = job_run_ids[0][0]
    results.append({
        "level": level, "source_run_id": run_id, "extension_run_id": new_run_id,
        "rpm": rpm, "theta": theta, "t_dump": t_dump, "walltime": walltime,
        "job_id": job_run_ids[0][1],
    })
    print(f"L{level} {run_id} -> extension {new_run_id} (job {job_run_ids[0][1]}, "
          f"rpm={rpm:g} theta={theta:g}, walltime={walltime})")

out = PROJECT_ROOT / "experiments" / "settling_extension_manifest.json"
out.write_text(json.dumps(results, indent=2))
print(f"\nSaved manifest: {out}")
