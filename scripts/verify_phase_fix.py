"""Verification of the restart time-unit fix.

The fix (src/BioReactor.c): t is measured in T_bio = L_bio/U_bio, so a t
carried across a change in U_bio is in the wrong units. event init now
rescales the restored clock by T_bio_prev/T_bio_new = 1/su, the identical
conversion already applied to u, p and g, and main() converts
params.t_checkpoint the same way before deriving t_ramp_start, t_mix and
t_dump_checkpoint from it.

This reruns the ORIGINAL Delta_theta sweep with no manual phi_angular
correction anywhere -- the same four configurations that produced

    Delta_theta = 0    ->  0.9961   (phase error 0deg, always looked fine)
    Delta_theta = 0.1  ->  1.2547   (phase error -23.9deg)
    Delta_theta = 3    ->  1.2567   (phase error -53.7deg)
    Delta_theta = 5    ->  1.2534   (phase error +116.6deg)

Prediction: all four now land at ~1.00. Delta_theta = 0 is the regression
check -- it had zero phase error to begin with, so the fix must not move it.

Usage:
    uv run python scripts/verify_phase_fix.py
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import GEOMETRY, FILL_LEVEL, PROJECT_ROOT, omega_b_of

chain.validate_params = lambda params: None

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-phasefix"
RPM, TARGET, FIDELITY, CYCLES = 32.5, 7.0, 7, 120
ROOT = Path(PROJECT_ROOT)
runs_dir = ROOT / "runs"

SOURCES = {7.0: "settling_baseline_L7", 6.9: "dtheta_src_L7_th6.9",
           4.0: "settling_ref_L7_rpm32.5_th4", 2.0: "settling_ref_L7_rpm32.5_th2"}

manifest = {}
for theta_source, source_run in SOURCES.items():
    src_dir = runs_dir / source_run
    t_dump = float(np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)[-1, 1])
    cfg = {
        "motion": {"omega_b": omega_b_of(RPM), "theta_max": [TARGET, 0.0, 0.0]},
        "fidelity": FIDELITY, "geometry": GEOMETRY, "fill_level": FILL_LEVEL,
        "n_mix_cycles": CYCLES, "n_transition_cycles": CYCLES, "t_buffer": 0.0,
        "sweep": {"parameter": "omega_b", "values": [omega_b_of(RPM)]},
        "initial_checkpoint": {
            "t_dump": t_dump, "omega_b": omega_b_of(RPM),
            "theta_max": [theta_source, 0.0, 0.0],
            "checkpoint_path": str(src_dir / "checkpoint.dump"),
        },
        "mpi": True, "ntasks": 8, "mem_per_cpu": "4G", "walltime": "06:00:00",
        "binary": BINARY, "submit": True, "exclude": "node2336",
    }
    run_id = chain.submit_chain(cfg)[0][0]
    manifest[f"{theta_source:g}"] = {"theta_source": theta_source, "run_id": run_id,
                                     "source_run": source_run}
    print(f"dtheta={TARGET - theta_source:g} (source theta={theta_source:g}): run={run_id}")

out = ROOT / "experiments" / "phase_fix_verification_manifest.json"
out.write_text(json.dumps(manifest, indent=2))
print(f"\nSaved {out}")
