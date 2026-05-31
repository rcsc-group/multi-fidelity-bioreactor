"""Submit restart_nc10: same as fix29_verify2 but n_mix_cycles=10.

T_per_st (omega_b=2pi, a=0.25, b=0.071) = 0.608 non-dim.
n_mix_cycles=10 → t_mix = t_ck + 6.08  (vs 1.84 for nc=3)
t_end_rel=10.0  → kLa window ≈ 10.0 - 6.08 = 3.92 non-dim  (vs 3.18 for nc=3)
"""
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
import simulate

CHECKPOINT = str(
    PROJECT_ROOT / "runs" / "5d4c73d7" / "checkpoint.dump"
)

params = {
    "run_id": "restart_nc10",
    "fidelity": 3,
    "geometry": {"a": 0.25, "b": 0.071, "n": 2.0},
    "fill_level": 0.5,
    "omega_b": 6.28318,
    "omega_b_prev": 3.14159,
    "theta_max": [7.0],
    "theta_max_prev": [7.0],
    "n_mix_cycles": 10,
    "t_checkpoint": 13.166515513678988,
    "t_end": 10.0,
    "restart_mode": "checkpoint",
}

runs_root = PROJECT_ROOT / "runs"
job_id = simulate.submit_slurm(
    params,
    project_root=PROJECT_ROOT,
    runs_root=runs_root,
    walltime="01:00:00",
    checkpoint=CHECKPOINT,
)
print(f"Submitted job {job_id}")
run_dir = runs_root / params["run_id"]
print(f"Run dir: {run_dir}")

results = simulate.wait_for_result(run_dir, timeout=1800, poll=15)
print(f"Results: {results}")
