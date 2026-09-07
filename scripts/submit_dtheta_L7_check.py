"""Resolution-independence check for the Delta_theta_max settling-time
pattern (diary.md 2026-09-07): does the near-flat/steep-onset behavior
found at L6 (target theta=7deg) replicate at L7, or does it change with
grid resolution (which would point toward a numerical/discretization
artifact rather than a real phenomenon)?

Reuses existing L7 cold-start references at no extra cost: settling_
baseline_L7 (theta=7, gives the trivial Delta_theta=0 case) and
settling_ref_L7_rpm32.5_th{2,4} (give Delta_theta=5 and 3 for free, from
the original 3x3 grid study). Only ONE new L7 source is needed: theta=6.9,
to test the critical small-Delta_theta=0.1 point where L6 showed the
sharpest part of the climb.

Same binary as the original L6 sweep (mainline, NOT the dAk/dt-fixed one)
for an apples-to-apples resolution comparison of the ORIGINAL phenomenon,
decoupled from the (already-shown-ineffective) ramp fix.

Usage:
    uv run python scripts/submit_dtheta_L7_check.py            # phase 1: new source
    uv run python scripts/submit_dtheta_L7_check.py transitions # phase 2, after phase 1 completes
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm
from scripts.settling_study_grid import T_END, GEOMETRY, FILL_LEVEL, PROJECT_ROOT, TEMPLATE, MAINLINE_BINARY, omega_b_of
import scripts.chain as chain

chain.validate_params = lambda params: None

RPM = 32.5
TARGET_THETA = 7.0
FIDELITY = 7
NTASKS = 8
CYCLES_PER_RUN = 60
THETA_SOURCES = [7.0, 6.9, 4.0, 2.0]  # d_theta = 0, 0.1, 3, 5


def source_run_id(theta: float) -> str:
    if theta == 7.0:
        return "settling_baseline_L7"
    if theta == 4.0:
        return "settling_ref_L7_rpm32.5_th4"
    if theta == 2.0:
        return "settling_ref_L7_rpm32.5_th2"
    return f"dtheta_src_L7_th{theta:g}"


def submit_source():
    theta = 6.9
    params = {
        "run_id": source_run_id(theta), "fidelity": FIDELITY, "geometry": GEOMETRY,
        "fill_level": FILL_LEVEL, "n_harmonics": 1, "theta_max": [theta, 0.0, 0.0],
        "phi_angular": [0.0, 0.0, 0.0], "omega_b": omega_b_of(RPM), "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": T_END, "n_mix_cycles": 80, "_binary": MAINLINE_BINARY,
    }
    job_id = submit_slurm(params, project_root=Path(PROJECT_ROOT), walltime="01:00:00",
                           template=Path(TEMPLATE), ntasks=NTASKS, mem="4G")
    print(f"L7 theta={theta:g} -> run={params['run_id']}  job={job_id}")


def submit_transitions():
    runs_dir = Path(PROJECT_ROOT) / "runs"
    manifest = {"target_theta": TARGET_THETA, "rpm": RPM, "fidelity": FIDELITY, "edges": {}}
    for theta_source in THETA_SOURCES:
        src_dir = runs_dir / source_run_id(theta_source)
        d = np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)
        t_dump = float(d[-1, 1])
        cfg = {
            "motion": {"omega_b": omega_b_of(RPM), "theta_max": [TARGET_THETA, 0.0, 0.0]},
            "fidelity": FIDELITY, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
            "n_mix_cycles": CYCLES_PER_RUN, "n_transition_cycles": CYCLES_PER_RUN, "t_buffer": 0.0,
            "sweep": {"parameter": "omega_b", "values": [omega_b_of(RPM)]},
            "initial_checkpoint": {
                "t_dump": t_dump, "omega_b": omega_b_of(RPM),
                "theta_max": [theta_source, 0.0, 0.0],
                "checkpoint_path": str(src_dir / "checkpoint.dump"),
            },
            "mpi": True, "ntasks": NTASKS, "mem_per_cpu": "4G", "walltime": "06:00:00",
            "binary": MAINLINE_BINARY, "submit": True,
        }
        job_run_ids = chain.submit_chain(cfg)
        run_id = job_run_ids[0][0]
        manifest["edges"][f"{theta_source:g}"] = {
            "theta_source": theta_source, "run_id": run_id, "source_run": source_run_id(theta_source),
        }
        print(f"theta_source={theta_source:g} -> target={TARGET_THETA:g}: run={run_id}  {job_run_ids}")

    out_path = Path(PROJECT_ROOT) / "experiments" / "dtheta_L7_check_manifest.json"
    out_path.write_text(json.dumps(manifest, indent=2))
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "transitions":
        submit_transitions()
    else:
        submit_source()
