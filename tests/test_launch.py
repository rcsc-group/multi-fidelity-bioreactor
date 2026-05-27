"""Tests for launch.py — SLURM script generation and run directory setup."""
import json
from pathlib import Path

import pytest

from scripts.launch import main as launch_main

VALID_PARAMS = {
    "fidelity": 7,
    "omega_b": 3.93,
    "n_harmonics": 2,
    "theta_max": [7.0, 2.0, 0.0],
    "phi_angular": [0.0, 1.57, 0.0],
    "omega_h": 2.1,
    "amplitude_h": [0.01, 0.005, 0.0],
    "phi_horizontal": [0.0, 0.78, 0.0],
    "geometry": {"a": 0.25, "b": 0.071, "n": 2.0},
    "fill_level": 0.5,
}


def test_slurm_script_contains_fidelity_level(tmp_path):
    """Generated SLURM script contains the LEVEL matching params fidelity."""
    params_file = tmp_path / "params.json"
    params_file.write_text(json.dumps(VALID_PARAMS))

    result = launch_main(str(params_file), runs_root=str(tmp_path / "runs"))

    slurm_script = Path(result["slurm_script"])
    assert f"LEVEL={VALID_PARAMS['fidelity']}" in slurm_script.read_text()


def test_run_directory_created(tmp_path):
    """launch.main() creates runs/{run_id}/ directory."""
    params_file = tmp_path / "params.json"
    params_file.write_text(json.dumps(VALID_PARAMS))

    result = launch_main(str(params_file), runs_root=str(tmp_path / "runs"))

    assert Path(result["run_dir"]).is_dir()


def test_params_json_copied_to_run_directory(tmp_path):
    """runs/{run_id}/params.json exists and matches the input params."""
    params_file = tmp_path / "params.json"
    params_file.write_text(json.dumps(VALID_PARAMS))

    result = launch_main(str(params_file), runs_root=str(tmp_path / "runs"))

    stored = json.loads((Path(result["run_dir"]) / "params.json").read_text())
    for key in VALID_PARAMS:
        assert stored[key] == VALID_PARAMS[key]


def test_unique_run_ids(tmp_path):
    """Two calls with the same params produce different run_ids."""
    params_file = tmp_path / "params.json"
    params_file.write_text(json.dumps(VALID_PARAMS))

    r1 = launch_main(str(params_file), runs_root=str(tmp_path / "runs"))
    r2 = launch_main(str(params_file), runs_root=str(tmp_path / "runs"))

    assert r1["run_id"] != r2["run_id"]
