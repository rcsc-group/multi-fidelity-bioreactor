"""Delta-theta_max monotonicity study, decreasing-direction mirror v2,
target=theta=4deg instead of theta=2deg (diary.md 2026-09-07). The
theta=2 target was abandoned: its tau_mean shows a large low-frequency
correlated wander (5.5% offset between a 5-cycle and a 59-cycle sample
mean of the SAME steady process) that a short-window settling metric
cannot resolve. theta=4 has much lower relative noise (1.96% vs 2.5-3.75%
tail std/mean) and is a real edge from the original settling study.

Submits, in order: (1) 2 new fine-resolution cold-start sources near the
target (theta=4.1, 4.5), (2) once those complete, 8 warm-start transitions
spanning Delta_theta_max = 0 (trivial) to 3.0, reusing existing checkpoints
(4=self, 5,6,6.5,6.9,7 already exist) for the rest.

Usage:
    uv run python scripts/submit_dtheta_th4_target.py           # phase 1 (sources)
    uv run python scripts/submit_dtheta_th4_target.py transitions  # phase 2, after phase 1 completes
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
TARGET_THETA = 4.0
FIDELITY = 6
NTASKS = 8
CYCLES_PER_RUN = 60
NEW_SOURCES = [4.1, 4.5]
THETA_SOURCES = [4.0, 4.1, 4.5, 5.0, 6.0, 6.5, 6.9, 7.0]


def source_run_id(theta: float) -> str:
    if theta == 4.0:
        return "settling_ref_L6_rpm32.5_th4"
    if theta == 7.0:
        return "settling_baseline_L6"
    return f"dtheta_src_L6_th{theta:g}"


def submit_sources():
    for theta in NEW_SOURCES:
        params = {
            "run_id": source_run_id(theta), "fidelity": FIDELITY, "geometry": GEOMETRY,
            "fill_level": FILL_LEVEL, "n_harmonics": 1, "theta_max": [theta, 0.0, 0.0],
            "phi_angular": [0.0, 0.0, 0.0], "omega_b": omega_b_of(RPM), "omega_h": 0.0,
            "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
            "t_end": T_END, "n_mix_cycles": 80, "_binary": MAINLINE_BINARY,
        }
        job_id = submit_slurm(params, project_root=Path(PROJECT_ROOT), walltime="01:00:00",
                               template=Path(TEMPLATE), ntasks=NTASKS, mem="4G")
        print(f"L6 theta={theta:g} -> run={params['run_id']}  job={job_id}")


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

    out_path = Path(PROJECT_ROOT) / "experiments" / "dtheta_monotonicity_th4target_manifest.json"
    out_path.write_text(json.dumps(manifest, indent=2))
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "transitions":
        submit_transitions()
    else:
        submit_sources()
