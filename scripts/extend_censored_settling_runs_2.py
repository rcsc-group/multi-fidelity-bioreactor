"""Second-round extension: of the 7 conditions extended by
scripts/extend_censored_settling_runs.py, only 2 showed genuine ongoing
decay past 120 cycles rather than a plateau matching their own reference's
intrinsic noise floor (see diary/conversation 2026-09-12: L6 17.5/th4, L8
37.5/th7, L7 17.5/th2 all plateaued AT their reference's own self-noise
level -- already converged, extending them further would be wasted
compute). L6 32.5/th4 and L7 32.5/th2 were still measurably above their
reference's noise floor at cycle 119, so extend those two by another 60
cycles (continuing from THEIR OWN checkpoint, i.e. a third same-condition
segment) to see whether they converge or plateau at a genuinely different
level.

Usage:
    uv run python scripts/extend_censored_settling_runs_2.py
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

# (level, run_id [2nd-round extension run], rpm, theta, walltime)
CANDIDATES = [
    (6, "fd00e7b8", 32.5, 4.0, "00:30:00"),   # 1st extension elapsed 00:02:24
    (7, "e2a0fe1e", 32.5, 2.0, "01:00:00"),   # 1st extension elapsed 00:08:42
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

    t_dump = p["t_checkpoint"] + p["t_end"]

    cfg = {
        "motion": {"omega_b": omega_b_of(rpm), "theta_max": [theta, 0.0, 0.0]},
        "fidelity": level,
        "geometry": GEOMETRY,
        "fill_level": FILL_LEVEL,
        "n_mix_cycles": EXTEND_CYCLES,
        "n_transition_cycles": EXTEND_CYCLES,
        "t_buffer": 0.0,
        "sweep": {"parameter": "omega_b", "values": [omega_b_of(rpm)]},
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

out = PROJECT_ROOT / "experiments" / "settling_extension_manifest_2.json"
out.write_text(json.dumps(results, indent=2))
print(f"\nSaved manifest: {out}")
