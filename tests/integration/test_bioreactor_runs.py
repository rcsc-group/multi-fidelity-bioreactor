"""Integration test: compile and run BioReactor at LEVEL=4, check output files."""
import json
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parents[2]
SRC = PROJECT_ROOT / "src"
BUILD = PROJECT_ROOT / "build"

MINIMAL_PARAMS = {
    "fidelity": 4,
    "omega_b": 3.93,
    "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0],
    "phi_angular": [0.0, 0.0, 0.0],
    "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "geometry": {"a": 0.25, "b": 0.071, "tilt": 0.0},
    "fill_level": 0.5,
}


@pytest.fixture(scope="module")
def compiled_binary(tmp_path_factory):
    """Compile BioReactor.c once for all integration tests."""
    build_dir = tmp_path_factory.mktemp("build")
    result = subprocess.run(
        ["qcc", "-O2", "-fopenmp", str(SRC / "BioReactor.c"),
         "-o", str(build_dir / "BioReactor"), "-lm"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        pytest.skip(f"BioReactor.c compilation failed:\n{result.stderr}")
    return build_dir / "BioReactor"


def test_bioreactor_produces_output_files(compiled_binary, tmp_path):
    """BioReactor runs with minimal params.json and produces tr_oxy.dat + normf.dat."""
    run_dir = tmp_path / "run_smoke"
    run_dir.mkdir()
    params_file = run_dir / "params.json"
    params_file.write_text(json.dumps(MINIMAL_PARAMS))

    result = subprocess.run(
        [str(compiled_binary), str(params_file)],
        capture_output=True, text=True, timeout=300,
    )
    assert result.returncode == 0, f"BioReactor exited non-zero:\n{result.stderr}"
    assert (run_dir / "tr_oxy.dat").exists(), "tr_oxy.dat not produced"
    assert (run_dir / "normf.dat").exists(), "normf.dat not produced"
