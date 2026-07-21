"""Runnable companion to the "Your first optimization loop" tutorial.

Runs the real multi-fidelity BO loop (real KRR-LR-GPR surrogate training,
real EI acquisition, real ExperimentData bookkeeping) against a synthetic
benchmark function instead of a real Basilisk simulation, so it finishes in
seconds with no build step and no SLURM access needed.

Usage:
    uv run python examples/synthetic_bo_demo.py
"""
from __future__ import annotations

import json
import math
import shutil
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.loop import run_loop

# The "true" optimum this demo is hunting for -- known in advance because we
# wrote the objective, which is what lets this tutorial show the loop
# actually converging rather than just "running without crashing".
OPT_OMEGA_B = 4.0
OPT_FILL = 0.5
LF_BIAS = -0.15  # low-fidelity systematically under-reports vs high-fidelity


def synthetic_kla(params: dict, lf_fidelity: int) -> float:
    """A smooth bump with a known peak; low-fidelity evaluations are biased low."""
    omega_b = params["omega_b"]
    fill = params["fill_level"]
    value = math.exp(-((omega_b - OPT_OMEGA_B) ** 2) / 2.0
                      - ((fill - OPT_FILL) ** 2) / 0.02)
    if params.get("fidelity") == lf_fidelity:
        value += LF_BIAS
    return value


def fake_submit_slurm(params, project_root=None, runs_root=None, **kwargs):
    run_dir = Path(runs_root) / params["run_id"]
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "params.json").write_text(json.dumps(params))
    return "demo-job"


def make_fake_wait_for_result(lf_fidelity: int):
    def fake_wait_for_result(run_dir, timeout=None, poll=None):
        params = json.loads((Path(run_dir) / "params.json").read_text())
        return {"kLa_25": synthetic_kla(params, lf_fidelity)}
    return fake_wait_for_result


def main() -> None:
    lf_fidelity, hf_fidelity = 5, 7
    with tempfile.TemporaryDirectory() as tmp:
        tmp_root = Path(tmp)
        shutil.copytree(PROJECT_ROOT / "config", tmp_root / "config")

        cfg = {
            "experiment_dir": "bo_demo",
            "lf_fidelity": lf_fidelity,
            "hf_fidelity": hf_fidelity,
            "n_lf_init": 4,
            "n_hf_init": 2,
            "n_iter": 4,
            "kla_key": "kLa_25",
            "n_candidates": 200,
            "walltime": "00:01:00",
            "job_timeout": 30,
        }

        with patch("scripts.simulate.submit_slurm", side_effect=fake_submit_slurm), \
             patch("scripts.simulate.wait_for_result",
                   side_effect=make_fake_wait_for_result(lf_fidelity)), \
             patch("scripts.simulate.get_job_elapsed_seconds", return_value=1.0), \
             patch("scripts.loop._PROJECT_ROOT", tmp_root):
            run_loop(cfg)

    print(f"\nTrue optimum was omega_b={OPT_OMEGA_B}, fill_level={OPT_FILL} "
          f"(kLa_25=1.0 at high fidelity).")


if __name__ == "__main__":
    main()
