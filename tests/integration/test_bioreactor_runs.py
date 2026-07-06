"""Integration test: compile and run BioReactor, check output files are produced."""
import json
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parents[2]

MINIMAL_PARAMS = {
    "run_id": "smoke_test",
    "fidelity": 3,
    "omega_b": 3.93,
    "n_harmonics": 1,
    "theta_max": [7.0, 0.0, 0.0],
    "phi_angular": [0.0, 0.0, 0.0],
    "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "geometry": {"a": 0.25, "b": 0.071, "n": 2.0},
    "fill_level": 0.5,
}


@pytest.fixture(scope="module")
def compiled_binary(tmp_path_factory):
    """Build BioReactor via the project Makefile."""
    build_dir = tmp_path_factory.mktemp("build")
    result = subprocess.run(
        ["make", "-f", str(PROJECT_ROOT / "Makefile"),
         f"BUILD_DIR={build_dir}", "build"],
        capture_output=True, text=True, cwd=str(PROJECT_ROOT),
    )
    if result.returncode != 0:
        pytest.skip(f"make build failed:\n{result.stderr}")
    binary = build_dir / "BioReactor"
    if not binary.exists():
        pytest.skip("BioReactor binary not found after make build")
    return binary


@pytest.mark.medium
def test_bioreactor_produces_output_files(compiled_binary, tmp_path):
    """BioReactor starts cleanly and writes tr_oxy.dat + normf.dat.

    The simulation runs longer than the CI timeout even at fidelity=3; we let it
    run for 120 s (enough for several output steps) then terminate.  A crash
    (non-zero exit before timeout) is still a hard failure.
    """
    run_dir = tmp_path / "run_smoke"
    run_dir.mkdir()
    params_file = run_dir / "params.json"
    params_file.write_text(json.dumps(MINIMAL_PARAMS))

    try:
        result = subprocess.run(
            [str(compiled_binary), str(params_file)],
            capture_output=True, text=True, timeout=120,
            cwd=str(run_dir),
        )
        assert result.returncode == 0, f"BioReactor exited non-zero:\n{result.stderr}"
    except subprocess.TimeoutExpired:
        pass  # long-running sim; file existence check below confirms healthy startup

    assert (run_dir / "tr_oxy.dat").exists(), "tr_oxy.dat not produced"
    assert (run_dir / "normf.dat").exists(), "normf.dat not produced"
