"""Falsification test for the cold-start-beats-warm-start anomaly
(diary.md 2026-09-07). Hypothesis under test: the slowness is
DIRECTIONAL -- warm-starting from a stronger-forcing condition (theta=7
or 4) into a weaker one (theta=2) is slow because the flow has to shed
excess kinetic energy/structure built up under the stronger forcing, a
process slower than a cold start's monotonic build-up in lockstep with
its own ramp. If true, the REVERSE direction -- warm-starting from
theta=2 (a state with genuinely LESS energy than the target) into a
stronger theta=4 or theta=7 -- should NOT show the same anomaly: it
should settle comparably to, or faster than, cold start, because there
is no excess to shed.

This directly falsifies (or supports) the hypothesis with zero
ambiguity: if reverse-direction warm start is also much slower than
cold start, "shedding excess energy" is wrong and something else (e.g.
a checkpoint-restore correctness bug) is more likely.

Reuses the existing theta=2, rpm=32.5 cold-start references
(settling_ref_L{6,7,8}_rpm32.5_th2) as checkpoint sources -- zero new
compute for the source runs, only the 6 new transition runs (2 edges x
3 levels) need to be submitted.

Usage:
    uv run python scripts/submit_settling_study_falsification.py
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import (
    LEVELS, PROJECT_ROOT, MAINLINE_BINARY, GEOMETRY, FILL_LEVEL, omega_b_of,
)

chain.validate_params = lambda params: None

RUNS_DIR = Path(PROJECT_ROOT) / "runs"
CYCLES_PER_RUN = 60  # same window as v2, for a like-for-like settling-time comparison

SOURCE = {"rpm": 32.5, "theta": 2.0}
SOURCE_RUN_ID = {6: "settling_ref_L6_rpm32.5_th2", 7: "settling_ref_L7_rpm32.5_th2", 8: "settling_ref_L8_rpm32.5_th2"}

EDGES = [
    (SOURCE, {"rpm": 32.5, "theta": 4.0}),  # reverse of the already-measured 4->2 edge
    (SOURCE, {"rpm": 32.5, "theta": 7.0}),  # reverse of the two-step 7->4->2 path, direct jump
]


def edge_tag(frm, to):
    return f"{frm['rpm']:g}_{frm['theta']:g}__{to['rpm']:g}_{to['theta']:g}"


manifest = {"cycles_per_run": CYCLES_PER_RUN, "edges": {}}

for fidelity, ntasks in LEVELS.items():
    src_dir = RUNS_DIR / SOURCE_RUN_ID[fidelity]
    d = np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)
    t_dump = float(d[-1, 1])

    for frm, to in EDGES:
        cfg = {
            "motion": {"omega_b": omega_b_of(to["rpm"]), "theta_max": [to["theta"], 0.0, 0.0]},
            "fidelity": fidelity,
            "geometry": GEOMETRY,
            "fill_level": FILL_LEVEL,
            "n_mix_cycles": CYCLES_PER_RUN,
            "n_transition_cycles": CYCLES_PER_RUN,
            "t_buffer": 0.0,
            "sweep": {"parameter": "omega_b", "values": [omega_b_of(to["rpm"])]},
            "initial_checkpoint": {
                "t_dump": t_dump,
                "omega_b": omega_b_of(frm["rpm"]),
                "theta_max": [frm["theta"], 0.0, 0.0],
                "checkpoint_path": str(src_dir / "checkpoint.dump"),
            },
            "mpi": True,
            "ntasks": ntasks,
            "mem_per_cpu": "4G",
            "walltime": "12:00:00",
            "binary": MAINLINE_BINARY,
            "submit": True,
        }
        job_run_ids = chain.submit_chain(cfg)
        run_id = job_run_ids[0][0]
        tag = edge_tag(frm, to)
        manifest["edges"].setdefault(str(fidelity), {})[tag] = {
            "from": frm, "to": to, "run_id": run_id, "source_run": SOURCE_RUN_ID[fidelity],
        }
        print(f"L{fidelity} {frm} -> {to}: {job_run_ids}")

out_path = Path(PROJECT_ROOT) / "experiments" / "settling_study_falsification_manifest.json"
out_path.write_text(json.dumps(manifest, indent=2))
print(f"\nSaved {out_path}")
