import sys
import json
import subprocess
import pathlib
import numpy as np
import pytest

# Make scripts/ importable from tests/
sys.path.insert(0, str(pathlib.Path(__file__).parents[1]))

PROJECT_ROOT = pathlib.Path(__file__).parents[1]
BINARY = PROJECT_ROOT / "build" / "BioReactor"

CANONICAL_PARAMS = {
    "run_id": "canonical",
    "fidelity": 3,
    "omega_b": 3.93,
    "n_harmonics": 1,
    "theta_max":      [7.0, 0.0, 0.0],
    "phi_angular":    [0.0, 0.0, 0.0],
    "omega_h": 0.0,
    "amplitude_h":    [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
    "fill_level": 0.5,
}


def pytest_configure(config):
    config.addinivalue_line("markers", "medium: short CFD run, fidelity 3; ~2 min on OSCAR")
    config.addinivalue_line("markers", "hpc: SLURM-required; fidelity 5-7")


def run_bioreactor(params: dict, tmp_path: pathlib.Path, timeout: int = 300) -> pathlib.Path:
    """Write params.json into a fresh run dir, execute BioReactor, return run_dir.

    Does not assert returncode — caller is responsible.  TimeoutExpired is caught
    silently so tests can check file existence regardless of sim completing.
    """
    if not BINARY.exists():
        pytest.skip(f"BioReactor binary not found at {BINARY}; run 'make build' first")
    run_dir = tmp_path / params["run_id"]
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "params.json").write_text(json.dumps(params))
    try:
        subprocess.run(
            [str(BINARY.resolve()), "params.json"],
            cwd=run_dir, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        pass
    return run_dir


def load_vol_frac(run_dir: pathlib.Path) -> np.ndarray:
    """Parse vol_frac_interf.dat → array (N, 6): i t f_liq_sum f_liq_interf posY_max posY_min"""
    path = run_dir / "vol_frac_interf.dat"
    lines = [l for l in path.read_text().splitlines()
             if l.strip() and not l.strip().startswith("i")]
    return np.array([[float(x) for x in l.split()] for l in lines])


def load_tr_oxy(run_dir: pathlib.Path) -> np.ndarray:
    """Parse tr_oxy.dat → array (N, 12): i t oxy_liq_sum oxy_liq_sum2 ..."""
    path = run_dir / "tr_oxy.dat"
    lines = [l for l in path.read_text().splitlines()
             if l.strip() and not l.strip().startswith("i")]
    return np.array([[float(x) for x in l.split()] for l in lines])


def load_normf(run_dir: pathlib.Path) -> np.ndarray:
    """Parse normf.dat → array (N, 14): i t Omega_avg Omega_rms ... ux ... uy ..."""
    path = run_dir / "normf.dat"
    lines = [l for l in path.read_text().splitlines()
             if l.strip() and not l.strip().startswith("i")]
    return np.array([[float(x) for x in l.split()] for l in lines])
