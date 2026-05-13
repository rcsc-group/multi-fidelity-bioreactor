"""Verify that DIAGNOSTICS=1 binary writes pressure_diag.dat with healthy residuals.

Requires build/BioReactor-health (compiled with DIAGNOSTICS=1 via 'make build-health').
The pressure Poisson residual from centered.h's multigrid solve must stay < 1e-4 for
a healthy run — values larger than that indicate the incompressibility constraint
is not being satisfied within the expected iteration budget.
"""
from __future__ import annotations

import sys
import json
import subprocess
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[2]))
from tests.conftest import CANONICAL_PARAMS

PROJECT_ROOT = Path(__file__).parents[2]
HEALTH_BINARY = PROJECT_ROOT / "build" / "BioReactor-health"


def _run_health(params: dict, tmp_path: Path, timeout: int = 300) -> Path:
    if not HEALTH_BINARY.exists():
        pytest.skip(f"BioReactor-health not found at {HEALTH_BINARY}; run 'make build-health' first")
    run_dir = tmp_path / params["run_id"]
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "params.json").write_text(json.dumps(params))
    try:
        subprocess.run(
            [str(HEALTH_BINARY.resolve()), "params.json"],
            cwd=run_dir, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        pass
    return run_dir


@pytest.mark.medium
def test_pressure_diag_file_written(tmp_path):
    """BioReactor-health must produce pressure_diag.dat when DIAGNOSTICS=1."""
    params = {**CANONICAL_PARAMS, "run_id": "pdiag_exists", "t_end": 1.0}
    run_dir = _run_health(params, tmp_path)
    assert (run_dir / "pressure_diag.dat").exists(), (
        "pressure_diag.dat not found — DIAGNOSTICS=1 binary may not be compiled"
    )


@pytest.mark.medium
def test_pressure_residuals_below_threshold(tmp_path):
    """Multigrid pressure residuals must stay < 1e-4 for a healthy run.

    Basilisk's centered.h Poisson solver (project()) typically converges to
    ~1e-6 to 1e-8 per step at fidelity=3. Residuals above 1e-4 indicate
    either insufficient multigrid iterations or solver divergence.
    """
    params = {**CANONICAL_PARAMS, "run_id": "pdiag_residuals", "t_end": 2.0}
    run_dir = _run_health(params, tmp_path)

    pdiag = run_dir / "pressure_diag.dat"
    if not pdiag.exists():
        pytest.skip("pressure_diag.dat not written — check DIAGNOSTICS=1 build")

    lines = [l for l in pdiag.read_text().splitlines()
             if l.strip() and not l.strip().startswith("i")]
    assert len(lines) >= 5, f"Too few rows in pressure_diag.dat: {len(lines)}"

    rows = np.array([[float(x) for x in l.split()] for l in lines])
    mgp_resa = rows[:, 2]   # col 2: mgp final residual
    resa_max = mgp_resa.max()

    assert resa_max < 1e-4, (
        f"Max pressure residual {resa_max:.2e} exceeds 1e-4 — Poisson solver not converging"
    )
