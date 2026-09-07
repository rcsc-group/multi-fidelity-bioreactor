"""Test C -- direct field comparison at the restore instant (diary.md
2026-09-07). Uses the new restart_diagnostic.txt instrumentation
(statsf() summary stats on u.x, u.y, p, f) written BOTH right before
dump() (source run's pre-dump state) and immediately after restore()
succeeds, BEFORE any ramp/forcing runs (restarted run's post-restore
state).

If restore() faithfully reconstructs exactly what dump() wrote, the
post_restore diagnostic from ANY restart off a given checkpoint --
trivial (Delta_theta=0) or the already-anomalous Delta_theta=0.05 -- should
match that checkpoint's own pre_dump diagnostic to near machine precision,
REGARDLESS of what target theta is requested (the diagnostic is written
before any forcing/ramp code runs). If they DO match: restore() itself is
clean, and the anomaly must originate in what happens AFTER restore
(forcing/physics evolution), not in field reconstruction. If they DON'T
match even for the trivial case: that is direct evidence of a restore()
bug independent of anything ramp-related.

Three short L6 runs, new instrumented binary (BioReactor-mpi-testc):
  1. source: cold start at theta=4.1 -> writes pre_dump diagnostic
  2. trivial restart: source -> target=4.1 (Delta_theta=0)
  3. anomalous restart: source -> target=4.05 (Delta_theta=0.05, already
     shown to plateau ~4.9% off in the L6 Test A)
Restarts are deliberately minimal (n_mix_cycles=1) -- only the instant
right after restore() matters here, not the subsequent evolution.

Usage:
    uv run python scripts/test_C_restore_fidelity.py            # phase 1: source
    uv run python scripts/test_C_restore_fidelity.py restarts    # phase 2, after phase 1 completes
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm
from scripts.settling_study_grid import T_END, GEOMETRY, FILL_LEVEL, PROJECT_ROOT, TEMPLATE, omega_b_of
import scripts.chain as chain

chain.validate_params = lambda params: None

TESTC_BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-testc"
RPM = 32.5
FIDELITY = 6
NTASKS = 8
SOURCE_RUN_ID = "testc_source_th4.1"


def submit_source():
    params = {
        "run_id": SOURCE_RUN_ID, "fidelity": FIDELITY, "geometry": GEOMETRY,
        "fill_level": FILL_LEVEL, "n_harmonics": 1, "theta_max": [4.1, 0.0, 0.0],
        "phi_angular": [0.0, 0.0, 0.0], "omega_b": omega_b_of(RPM), "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": T_END, "n_mix_cycles": 80, "_binary": TESTC_BINARY,
    }
    job_id = submit_slurm(params, project_root=Path(PROJECT_ROOT), walltime="01:00:00",
                           template=Path(TEMPLATE), ntasks=NTASKS, mem="4G")
    print(f"source theta=4.1 -> run={SOURCE_RUN_ID}  job={job_id}")


def submit_restarts():
    runs_dir = Path(PROJECT_ROOT) / "runs"
    src_dir = runs_dir / SOURCE_RUN_ID
    d = np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)
    t_dump = float(d[-1, 1])

    for tag, target_theta in [("trivial", 4.1), ("anomalous", 4.05)]:
        cfg = {
            "motion": {"omega_b": omega_b_of(RPM), "theta_max": [target_theta, 0.0, 0.0]},
            "fidelity": FIDELITY, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
            "n_mix_cycles": 1, "n_transition_cycles": 1, "t_buffer": 0.0,
            "sweep": {"parameter": "omega_b", "values": [omega_b_of(RPM)]},
            "initial_checkpoint": {
                "t_dump": t_dump, "omega_b": omega_b_of(RPM),
                "theta_max": [4.1, 0.0, 0.0],
                "checkpoint_path": str(src_dir / "checkpoint.dump"),
            },
            "mpi": True, "ntasks": NTASKS, "mem_per_cpu": "4G", "walltime": "01:00:00",
            "binary": TESTC_BINARY, "submit": True,
        }
        job_run_ids = chain.submit_chain(cfg)
        print(f"{tag} (target={target_theta:g}): run={job_run_ids[0][0]}  {job_run_ids}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "restarts":
        submit_restarts()
    else:
        submit_source()
