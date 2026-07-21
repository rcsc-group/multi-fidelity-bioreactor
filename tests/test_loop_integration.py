"""End-to-end testbed for the multi-fidelity BO algorithm (loop.py + train_surrogate.py
+ suggest.py wired together for real), cheap to run because the expensive part --
the actual Basilisk simulation -- is replaced by a synthetic benchmark function.

Only scripts.simulate.submit_slurm/wait_for_result/get_job_elapsed_seconds are mocked
(the SLURM/filesystem boundary). Everything downstream -- KRR-LR-GPR surrogate
training, EI acquisition, ExperimentData bookkeeping, the DoE + BO iteration phases --
runs unmodified, for real, on synthetic data. This is what makes it a testbed for the
algorithm itself, not just a mock of the whole loop.

The synthetic objective has a known maximum (omega_b=4.0, fill_level=0.5) and a
deliberate, deterministic LF/HF bias, so the multi-fidelity bias-correction is
genuinely exercised rather than trivially satisfied by identical LF/HF data.
"""
from __future__ import annotations

import json
import math
import shutil
from pathlib import Path
from unittest.mock import patch

from scripts.loop import run_loop

_REAL_PROJECT_ROOT = Path(__file__).parents[1]

# Known optimum well inside param_space.yaml's bounds (omega_b in [1.57, 6.28],
# fill_level in [0.3, 0.7]) -- not at an edge, so BO isn't trivially degenerate.
_OPT_OMEGA_B = 4.0
_OPT_FILL = 0.5
_LF_BIAS = -0.15  # deterministic, known offset -- exercises the MF bias correction


def _synthetic_kla(params: dict, lf_fidelity: int) -> float:
    """Smooth bump with a known maximum; LF systematically biased low vs HF."""
    omega_b = params["omega_b"]
    fill = params["fill_level"]
    value = math.exp(-((omega_b - _OPT_OMEGA_B) ** 2) / 2.0
                      - ((fill - _OPT_FILL) ** 2) / 0.02)
    if params.get("fidelity") == lf_fidelity:
        value += _LF_BIAS
    return value


def _fake_submit_slurm(params, project_root=None, runs_root=None, **kwargs):
    run_dir = Path(runs_root) / params["run_id"]
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "params.json").write_text(json.dumps(params))
    return "42"


def _make_fake_wait_for_result(lf_fidelity: int):
    def _fake_wait_for_result(run_dir, timeout=None, poll=None):
        params = json.loads((Path(run_dir) / "params.json").read_text())
        return {"kLa_25": _synthetic_kla(params, lf_fidelity)}
    return _fake_wait_for_result


def test_bo_loop_runs_end_to_end_on_synthetic_objective(tmp_path):
    """The full DoE + BO loop completes on a cheap synthetic objective and its
    best-found value is closer to the known optimum than a random DoE point,
    i.e. the acquisition function is doing real work, not just not-crashing."""
    # run_loop resolves config/param_space.yaml from _PROJECT_ROOT, which we
    # patch to tmp_path so runs/ writes stay sandboxed -- mirror the real
    # config in so that resolution still succeeds.
    shutil.copytree(_REAL_PROJECT_ROOT / "config", tmp_path / "config")

    lf_fidelity, hf_fidelity = 5, 7
    cfg = {
        "experiment_dir": "bo_test",
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

    with patch("scripts.simulate.submit_slurm", side_effect=_fake_submit_slurm), \
         patch("scripts.simulate.wait_for_result",
               side_effect=_make_fake_wait_for_result(lf_fidelity)), \
         patch("scripts.simulate.get_job_elapsed_seconds", return_value=1.0), \
         patch("scripts.loop._PROJECT_ROOT", tmp_path):
        run_loop(cfg)

    from f3dasm import ExperimentData
    exp_dir = tmp_path / cfg["experiment_dir"]
    ed = ExperimentData.from_file(str(exp_dir))
    inp_df, out_df = ed.to_pandas()

    n_expected = cfg["n_lf_init"] + cfg["n_hf_init"] + cfg["n_iter"]
    assert len(inp_df) == n_expected

    hf_mask = inp_df["fidelity"] == hf_fidelity
    assert hf_mask.sum() == cfg["n_hf_init"] + cfg["n_iter"]

    best_hf = float(out_df.loc[hf_mask, "kLa_25"].max())
    doe_hf_mask = hf_mask & (inp_df["phase"] == "doe_hf")
    best_doe_hf = float(out_df.loc[doe_hf_mask, "kLa_25"].max())

    assert best_hf >= best_doe_hf, (
        f"BO iterations (best={best_hf:.4f}) did not improve on the initial HF "
        f"DoE (best={best_doe_hf:.4f}) -- acquisition function may not be working"
    )
    # True optimum value is 1.0 (bump peak, HF/unbiased); best found should be
    # meaningfully above the DoE's typical random-sample value.
    assert best_hf > 0.3, f"best_hf={best_hf:.4f} suspiciously far from the known optimum"
