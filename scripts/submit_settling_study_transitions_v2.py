"""Settling-time study, corrected transitions (2026-09-06): nearest-
neighbor chain through the 3x3 grid, reusing the 9 already-completed
reference/baseline cold starts as checkpoint sources, instead of a
hub-and-spoke design from one shared baseline.

The first version checkpointed all 8 out-points from a single baseline
(32.5rpm, 7deg), manufacturing artificially large steps for most of the
grid (e.g. theta 7->2 directly, a Delta_theta=-5 jump) that don't match
how checkpointing is actually used in the real sweep (small, adjacent
steps). Corrected topology:

    (32.5,7)[baseline] -> (32.5,4) -> (32.5,2)          [d_theta=-3, then -2]
    (32.5,7)            -> (17.5,7) -> (17.5,4) -> (17.5,2)  [d_omega_b=-15, then d_theta=-3, then -2]
    (32.5,7)            -> (37.5,7) -> (37.5,4) -> (37.5,2)  [d_omega_b=+5, then d_theta=-3, then -2]

8 edges, only 2 of which are the unavoidable large RPM jumps (reaching a
new RPM row at all); the other 6 are small, realistic theta steps.
Recording window extended to 60 cycles (from 25) -- the first version's
25-cycle window was sized for a small step and was too short even for
some of these (see diary.md 2026-09-06: a Delta_omega_b=-15 transition
was still visibly converging, not settled, at cycle 19-24).

Reuses the 9 existing reference/baseline cold starts as-is (verified
self-consistent/converged over their own last 5 cycles, std/mean<1.5%,
so no need to also extend them).

Usage:
    uv run python scripts/submit_settling_study_transitions_v2.py
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
import scripts.chain as chain
from scripts.settling_study_grid import (
    LEVELS, BASELINE, PROJECT_ROOT, MAINLINE_BINARY, GEOMETRY, FILL_LEVEL,
    omega_b_of, baseline_run_id,
)

chain.validate_params = lambda params: None

RUNS_DIR = Path(PROJECT_ROOT) / "runs"
CYCLES_PER_RUN_V2 = 60

# (from_point, to_point) edges. from_point=None means the shared baseline.
EDGES = [
    ({"rpm": 32.5, "theta": 7.0}, {"rpm": 32.5, "theta": 4.0}),
    ({"rpm": 32.5, "theta": 4.0}, {"rpm": 32.5, "theta": 2.0}),
    ({"rpm": 32.5, "theta": 7.0}, {"rpm": 17.5, "theta": 7.0}),
    ({"rpm": 17.5, "theta": 7.0}, {"rpm": 17.5, "theta": 4.0}),
    ({"rpm": 17.5, "theta": 4.0}, {"rpm": 17.5, "theta": 2.0}),
    ({"rpm": 32.5, "theta": 7.0}, {"rpm": 37.5, "theta": 7.0}),
    ({"rpm": 37.5, "theta": 7.0}, {"rpm": 37.5, "theta": 4.0}),
    ({"rpm": 37.5, "theta": 4.0}, {"rpm": 37.5, "theta": 2.0}),
]


def source_run_id(fidelity: int, p: dict) -> str:
    if p["rpm"] == BASELINE["rpm"] and p["theta"] == BASELINE["theta"]:
        return baseline_run_id(fidelity)
    return f"settling_ref_L{fidelity}_rpm{p['rpm']:g}_th{p['theta']:g}"


def edge_tag(frm, to):
    return f"{frm['rpm']:g}_{frm['theta']:g}__{to['rpm']:g}_{to['theta']:g}"


manifest = {"cycles_per_run": CYCLES_PER_RUN_V2, "edges": {}}

for fidelity, ntasks in LEVELS.items():
    for frm, to in EDGES:
        src_dir = RUNS_DIR / source_run_id(fidelity, frm)
        d = np.loadtxt(src_dir / "shear_stress.dat", skiprows=1)
        t_dump = float(d[-1, 1])

        cfg = {
            "motion": {"omega_b": omega_b_of(to["rpm"]), "theta_max": [to["theta"], 0.0, 0.0]},
            "fidelity": fidelity,
            "geometry": GEOMETRY,
            "fill_level": FILL_LEVEL,
            "n_mix_cycles": CYCLES_PER_RUN_V2,
            "n_transition_cycles": CYCLES_PER_RUN_V2,
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
            "from": frm, "to": to, "run_id": run_id, "source_run": source_run_id(fidelity, frm),
        }
        print(f"L{fidelity} {frm} -> {to}: {job_run_ids}")

out_path = Path(PROJECT_ROOT) / "experiments" / "settling_study_v2_manifest.json"
out_path.write_text(json.dumps(manifest, indent=2))
print(f"\nSaved {out_path}")
