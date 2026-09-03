"""Does the same-condition restart-chain bias result (diary.md 2026-09-03
(2): tau_mean_max +0.1%, tau_100_max -4.0% at 32.5rpm) hold across other
conditions, or was it specific to that one RPM?

Same two-arm design as submit_restart_bias_test_l6.py (fresh cold start
vs. 4x6-cycle same-condition checkpoint chain, L6, BioReactor-mpi-
rampmatch), repeated at a spread of RPMs spanning the fig13a sweep:
17.5 (low), 25.0 (mid), 37.5 (high -- the condition with the worst
agreement vs. Kim and the longest ramp in cycles, so plausibly the
likeliest to show a bigger restart effect if one exists). 32.5rpm is not
resubmitted -- already have that data (runs/restart_bias_fresh_l6 +
459e1b76/d74d35b7/a3070b72/7beefe12).

Writes a manifest (runs/_restart_bias_multi_manifest.json) mapping each
RPM to its fresh run_id and ordered chain segment run_ids, since
chain.py's build_chain() assigns run_ids via uuid4 -- not predictable,
and calling build_chain() again separately would reseed them.

Usage:
    uv run python scripts/submit_restart_bias_multi.py
"""
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm, _t_mix_nd
import scripts.chain as chain

# Same validate_params bypass as submit_restart_bias_test_l6.py -- Kim's
# literal geometry sits outside config/param_space.yaml's BO search bounds
# by design; scoped to this script only.
chain.validate_params = lambda params: None

PROJECT_ROOT = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor"
TEMPLATE = f"{PROJECT_ROOT}/config/slurm_mpi_template.sh"
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-rampmatch"
MANIFEST_PATH = Path(PROJECT_ROOT) / "runs" / "_restart_bias_multi_manifest.json"

THETA_MAX = [7.0, 0.0, 0.0]
FIDELITY = 6
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
FILL_LEVEL = 0.5

N_SEGMENTS = 4
CYCLES_PER_SEGMENT = 6
FRESH_CYCLES_NOMINAL = 30

RPMS = [17.5, 25.0, 37.5]

manifest: dict[str, dict] = {}

for rpm in RPMS:
    omega_b = rpm * 2 * math.pi / 60.0
    base_params = {"omega_b": omega_b, "theta_max": THETA_MAX, "geometry": GEOMETRY, "n_mix_cycles": 1}
    T_per_nd = _t_mix_nd(base_params)
    print(f"--- {rpm}rpm: T_per_nd={T_per_nd:.6f} ---")

    fresh_run_id = f"restart_bias_fresh_l6_rpm{rpm:g}"
    fresh_params = {
        "run_id": fresh_run_id,
        "fidelity": FIDELITY,
        "geometry": GEOMETRY,
        "fill_level": FILL_LEVEL,
        "n_harmonics": 1,
        "theta_max": THETA_MAX,
        "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": omega_b,
        "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0],
        "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": FRESH_CYCLES_NOMINAL * T_per_nd,
        "n_mix_cycles": 80,
        "_binary": BINARY,
    }
    fresh_job = submit_slurm(
        fresh_params, project_root=PROJECT_ROOT, walltime="00:30:00",
        template=TEMPLATE, mem="2G", cpus=1, ntasks=4,
    )
    print(f"  Arm A (fresh): run_id={fresh_run_id} job={fresh_job}")

    chain_cfg = {
        "motion": {"omega_b": omega_b, "theta_max": THETA_MAX},
        "fidelity": FIDELITY,
        "geometry": GEOMETRY,
        "fill_level": FILL_LEVEL,
        "n_mix_cycles": CYCLES_PER_SEGMENT,
        "n_transition_cycles": CYCLES_PER_SEGMENT,
        "t_buffer": 0.0,
        "sweep": {"parameter": "omega_b", "values": [omega_b] * N_SEGMENTS},
        "mpi": True,
        "ntasks": 4,
        "mem_per_cpu": "2G",
        "walltime": "00:30:00",
        "binary": BINARY,
        "submit": True,
    }
    chain_run_job_ids = chain.submit_chain(chain_cfg)
    chain_run_ids = [run_id for run_id, _ in chain_run_job_ids]
    print(f"  Arm B (chained): {chain_run_job_ids}")

    manifest[f"{rpm:g}"] = {"fresh_run_id": fresh_run_id, "chain_run_ids": chain_run_ids}

MANIFEST_PATH.write_text(json.dumps(manifest, indent=2))
print(f"\nManifest written to {MANIFEST_PATH}")
