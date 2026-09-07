"""Re-run of the L7 Delta_theta_max check (diary.md 2026-09-07), now using
the p.nodump-restore-fixed binary (also has the dAk/dt ramp fix), to see
whether fixing the pressure-restore bug changes the pattern found with the
broken binary (target=theta=7deg, all 3 nonzero Delta_theta points never
settled in 60 cycles, plateauing at 12-20% offsets).

Reuses the EXACT SAME checkpoints as the original check (checkpoint.dump
content is unaffected by a binary-side fix -- only what happens AFTER
restore changes): settling_baseline_L7, dtheta_src_L7_th6.9,
settling_ref_L7_rpm32.5_th{4,2}. Window extended 60->120 cycles: the L6
verification of this same fix needed the full 60-cycle window just to
barely cross into tolerance (settling at cycle 57) for a much smaller
Delta_theta=0.05 case, so L7's larger offsets (12-20% pre-fix) plausibly
need more room to fully converge.

Usage:
    uv run python scripts/submit_dtheta_L7_check_v2_fixed.py
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

FIXED_BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-testc"
RPM = 32.5
TARGET_THETA = 7.0
FIDELITY = 7
NTASKS = 8
CYCLES_PER_RUN = 120  # doubled from 60 -- see docstring
THETA_SOURCES = [7.0, 6.9, 4.0, 2.0]  # d_theta = 0, 0.1, 3, 5


def source_run_id(theta: float) -> str:
    if theta == 7.0:
        return "settling_baseline_L7"
    if theta == 4.0:
        return "settling_ref_L7_rpm32.5_th4"
    if theta == 2.0:
        return "settling_ref_L7_rpm32.5_th2"
    return f"dtheta_src_L7_th{theta:g}"


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
        "binary": FIXED_BINARY, "submit": True,
    }
    job_run_ids = chain.submit_chain(cfg)
    run_id = job_run_ids[0][0]
    manifest["edges"][f"{theta_source:g}"] = {
        "theta_source": theta_source, "run_id": run_id, "source_run": source_run_id(theta_source),
    }
    print(f"theta_source={theta_source:g} -> target={TARGET_THETA:g}: run={run_id}  {job_run_ids}")

out_path = Path(PROJECT_ROOT) / "experiments" / "dtheta_L7_check_v2_fixed_manifest.json"
out_path.write_text(json.dumps(manifest, indent=2))
print(f"\nSaved {out_path}")
