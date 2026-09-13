"""Fig 8 replica investigation: our replica (t/Tp~34-37) is compared at a
much earlier point than Kim's own published window (t/Tp~80-83) -- never
tested as a hypothesis for the standing ~3-4x tau/EDR amplitude gap.
Extends l10_kim_fig8_signed's own checkpoint (t=22.46, t/Tp~37.0) by 25
more cycles (t/Tp~62.0, about 2/3 of the way to Kim's window) as a first,
bounded look at whether the amplitude is still trending toward Kim's
values -- not yet committing to the full ~46 cycles needed to reach t/Tp=83.

Same binary as the original run (BioReactor-mpi-video, still present in
scratch) -- avoids any new cross-binary risk for this specific check.
Trivial (same-condition, 32.5rpm/theta=7 unchanged) restart, so this
session's phase-offset-restart fix (only matters for a CONDITION change)
is moot here either way.

Walltime sized from this run's OWN first extension segment (job 4812614,
sacct Elapsed 01:27:49 for 2.98 cycles at 64 ranks = 29.5 min/cycle,
recorded in scripts/cost_model.py) x 1.8 safety margin -- one historical
data point at this exact (level, ntasks), so a generous margin, not the
minimum.

Usage:
    uv run python scripts/extend_fig8_late_time.py
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_model import t_per_nd

chain.validate_params = lambda params: None

PROJECT_ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
RUNS = PROJECT_ROOT / "runs"

SOURCE_RUN = "l10_kim_fig8_signed"
RPM, THETA = 32.5, 7.0
# [2026-09-12] Reduced from 25 to 10 cycles AND switched 64->48 ranks after
# job 6311172 sat PENDING with a predicted 2-day start: it requested 64
# tasks on a SINGLE node, a pool of exactly 4 batch-partition nodes, all 4
# unavailable (3 alloc, 1 drain) at submission time. 48-core nodes are a
# much larger pool (152 total, 61 with some free capacity). Slower per
# cycle (83.0 vs 29.5 min/cycle, scripts/cost_model.py) but should actually
# start soon instead of queuing behind a resource class we can't get.
# Shrunk the cycle count to keep total walltime reasonable against a hard
# 4-day deadline.
EXTEND_CYCLES = 10
# [2026-09-12, 2nd resubmit] 48 ranks got the SAME predicted 2-day start as
# 64 -- checked real-time free capacity across the whole batch partition:
# zero nodes have 48 idle cores right now (max free on any node was 33,
# then 38 a few minutes later). Genuine cluster-wide crunch, not a
# resource-shape mismatch. Shrunk to 32 ranks to fit an actual currently-
# free fragment (node2410 had 38 idle at last check). NO measured cost_model
# entry exists at (10, 32) -- walltime below is an EXPLICIT ESTIMATE
# (naive 1/ranks scaling off the 48-rank number, which was itself measured
# under contention so may already be pessimistic), not a measurement.
# Plan: check real throughput from logstats.dat within ~2h of it actually
# starting and be ready to react, rather than trusting this blind for the
# full walltime.
NTASKS = 32
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-video"

p = json.loads((RUNS / SOURCE_RUN / "params.json").read_text())
d_last = None
with open(RUNS / SOURCE_RUN / "shear_stress.dat") as fh:
    for line in fh:
        d_last = line
t_now = float(d_last.split()[1])
T = t_per_nd(RPM, THETA)
print(f"source run's last recorded t={t_now:.4f} (t/Tp={t_now/T:.2f})")

cfg = {
    "motion": {"omega_b": p["omega_b"], "theta_max": [THETA, 0.0, 0.0]},
    "fidelity": 10,
    "geometry": p["geometry"],
    "fill_level": p["fill_level"],
    "n_mix_cycles": EXTEND_CYCLES,
    "n_transition_cycles": EXTEND_CYCLES,
    "t_buffer": 0.0,
    "sweep": {"parameter": "omega_b", "values": [p["omega_b"]]},  # trivial: no change
    "initial_checkpoint": {
        "t_dump": t_now,
        "omega_b": p["omega_b"],
        "theta_max": [THETA, 0.0, 0.0],
        "checkpoint_path": str((RUNS / SOURCE_RUN / "checkpoint.dump").resolve()),
    },
    "mpi": True,
    "ntasks": NTASKS,
    "mem_per_cpu": "2G",
    "walltime": "36:00:00",  # ESTIMATE (no measurement at 32 ranks): 10 cycles x
                              # ~124.5 min/cycle (83.0 min/cyc @ 48 ranks x 48/32
                              # naive scaling) x ~1.75x extra margin for the
                              # scaling assumption itself being unverified
    "binary": BINARY,
    "videos": False,   # domain-mean time series only for this probe; re-run
                        # with field export once we know amplitude actually
                        # trends toward Kim's values
    "submit": True,
}
job_run_ids = chain.submit_chain(cfg)
new_run_id, job_id = job_run_ids[0]
print(f"submitted extension {new_run_id} (job {job_id}), target t/Tp="
      f"{(t_now + EXTEND_CYCLES * T)/T:.2f}")

out = PROJECT_ROOT / "experiments" / "fig8_late_time_extension_manifest.json"
out.write_text(json.dumps({
    "source_run_id": SOURCE_RUN, "extension_run_id": new_run_id, "job_id": job_id,
    "t_checkpoint": t_now, "extend_cycles": EXTEND_CYCLES,
}, indent=2))
print(f"saved {out}")
