"""Does same-condition checkpoint-restart chaining bias the result vs. a
single continuous cold start? (diary.md 2026-09-03)

This is the question the L10 fresh_mpi-vs-chain_mpi comparison was set up
to answer earlier in this project but never concluded (the two runs
ended up at non-overlapping cycle counts). It's also the open question
behind `runs/l10_kim_fig8_signed` (Fig 8a-c's data source): that run is
itself a same-condition restart continuation, chained through
l10_kim_seg0 -> seg1 -> seg2 -> fig8_signed, and we have no clean
cold-start baseline reaching that far to check it against.

Testing this at L10 is expensive (days/condition, needs multi-day
checkpoint chaining just to run one arm). Testing it at L6 is cheap and
exercises the exact same checkpoint-write/restart-read code path -- the
mechanism under suspicion doesn't care about grid resolution.

Two arms, same condition (32.5rpm/theta=7deg, matching Fig 8's own
condition and Kim's baseline), same total elapsed time (24 rocking
cycles), same driver (BioReactor-mpi-rampmatch, current bug-fixed
source):
  A. FRESH: one continuous cold start, 0 -> 24 cycles, no restarts.
  B. CHAINED: 4 segments of 6 cycles each, same omega_b/theta_max
     throughout (a "sweep" over 4 identical values -- chain.py doesn't
     care that the value doesn't change, it just chains via checkpoint
     restart 3 times), mimicking the seg0/seg1/seg2/fig8_signed hop
     count that produced the real Fig 8 data.

If A and B agree (within the ~1% noise floor already established for
tau_100/tau_mean via the OpenMP-race-fix reproducibility test), same-
condition restart chaining is ruled out as a bias source and
l10_kim_fig8_signed's data can be trusted as-is (just needs to be
extended further, at L10, for the real figure). If they diverge
meaningfully, that's a real, previously-uncharacterized bug in the
checkpoint-restart path, independent of resolution.

Usage:
    uv run python scripts/submit_restart_bias_test_l6.py
"""
import math
import os
import sys

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm, _t_mix_nd
import scripts.chain as chain

# chain.py's submit_chain() calls validate_params(), which enforces
# config/param_space.yaml -- the project's OWN Bayesian-optimization search
# domain (geometry.b in [0.05, 0.15]). Kim's literal geometry (b=0.03575)
# is deliberately outside that domain -- it's a fixed validation target,
# not a BO search point -- which is why every other Kim-replication script
# this session (fig13a_rampmatch, fig13a_l6) calls submit_slurm() directly
# and never goes through validate_params() at all. Same bypass here, scoped
# to this script only -- chain.py itself is left untouched so its real
# safety net for actual optimization sweeps stays intact.
chain.validate_params = lambda params: None
submit_chain = chain.submit_chain

PROJECT_ROOT = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor"
TEMPLATE = f"{PROJECT_ROOT}/config/slurm_mpi_template.sh"
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-rampmatch"

OMEGA_B = 32.5 * 2 * math.pi / 60.0   # Kim's baseline condition, matches Fig 8
THETA_MAX = [7.0, 0.0, 0.0]
FIDELITY = 6
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
FILL_LEVEL = 0.5

N_SEGMENTS = 4
CYCLES_PER_SEGMENT = 6
# NOTE: both chain.py and the C driver round each segment's actual stop time
# UP to the next full-period boundary (int(t_end_abs/T_per_st)+1) -- verified
# via a dry run that this adds +1 real cycle at EVERY segment boundary here
# (floating-point noise pushes the ratio just past the integer, not under
# it), so 4 nominal 6-cycle segments land at ~28 real cycles, not 24. Rather
# than predict the exact landing bit-for-bit, request the fresh arm well
# past it (30 nominal cycles, itself rounded up by the same C-side logic)
# and compare over the overlapping tail window using actual recorded t,
# not assumed cycle counts -- the same non-overlapping-window mistake that
# stalled the earlier fresh_mpi-vs-chain_mpi (L10) comparison.
FRESH_CYCLES_NOMINAL = 30

base_params = {
    "omega_b": OMEGA_B,
    "theta_max": THETA_MAX,
    "geometry": GEOMETRY,
    "n_mix_cycles": 1,
}
T_per_nd = _t_mix_nd(base_params)
print(f"T_per_nd = {T_per_nd:.6f}, fresh nominal t_end = {FRESH_CYCLES_NOMINAL * T_per_nd:.4f} "
      f"(will itself round up to the next period boundary)")

# ── Arm A: fresh, continuous, cold start ────────────────────────────────────
fresh_params = {
    "run_id": "restart_bias_fresh_l6",
    "fidelity": FIDELITY,
    "geometry": GEOMETRY,
    "fill_level": FILL_LEVEL,
    "n_harmonics": 1,
    "theta_max": THETA_MAX,
    "phi_angular": [0.0, 0.0, 0.0],
    "omega_b": OMEGA_B,
    "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "t_end": FRESH_CYCLES_NOMINAL * T_per_nd,
    "n_mix_cycles": 80,
    "_binary": BINARY,
}
DRY_RUN = bool(os.environ.get("DRY_RUN"))

chain_cfg = {
    "motion": {"omega_b": OMEGA_B, "theta_max": THETA_MAX},
    "fidelity": FIDELITY,
    "geometry": GEOMETRY,
    "fill_level": FILL_LEVEL,
    "n_mix_cycles": CYCLES_PER_SEGMENT,
    "n_transition_cycles": CYCLES_PER_SEGMENT,
    "t_buffer": 0.0,
    "sweep": {"parameter": "omega_b", "values": [OMEGA_B] * N_SEGMENTS},
    "mpi": True,
    "ntasks": 4,
    "mem_per_cpu": "2G",
    "walltime": "00:30:00",
    "binary": BINARY,
    "submit": not DRY_RUN,
}

if DRY_RUN:
    from scripts.chain import build_chain
    for k, p in enumerate(build_chain(chain_cfg)):
        print(f"seg {k}: run_id={p['run_id']} t_end={p['t_end']:.4f} "
              f"t_checkpoint={p.get('t_checkpoint')} omega_b_prev={p.get('omega_b_prev')} "
              f"_binary={p.get('_binary')}")
    print(f"fresh: t_end={fresh_params['t_end']:.4f} _binary={fresh_params['_binary']}")
    sys.exit(0)

if os.environ.get("SKIP_FRESH"):
    print("Arm A: skipped (already submitted -- job 5700520)")
else:
    fresh_job = submit_slurm(
        fresh_params,
        project_root=PROJECT_ROOT,
        walltime="00:30:00",
        template=TEMPLATE,
        mem="2G",
        cpus=1,
        ntasks=4,
    )
    print(f"Arm A (fresh, {FRESH_CYCLES_NOMINAL} nominal cycles): job {fresh_job}")

# ── Arm B: chained, 4 same-condition restart segments ───────────────────────
chain_job_ids = submit_chain(chain_cfg)
print(f"Arm B (chained, {N_SEGMENTS}x{CYCLES_PER_SEGMENT} cycles): jobs {chain_job_ids}")
