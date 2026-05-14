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


_NITERMAX = 1000  # must match BioReactor.c: NITERMAX = 1000


@pytest.mark.medium
def test_poisson_solver_converges(tmp_path):
    """Multigrid Poisson solver must not exhaust its iteration budget.

    project() uses tolerance=TOLERANCE/sq(dt), so mgp.resa is not comparable
    to a fixed threshold.  The correct health signal is mgp.i < NITERMAX=1000:
    hitting the budget means the divergence-free constraint was not satisfied
    and Basilisk would have printed a convergence warning to stderr.
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
    mgp_i = rows[:, 4]   # col 4: mgp_i (number of V-cycles)
    max_i = int(mgp_i.max())

    assert max_i < _NITERMAX, (
        f"Poisson solver hit iteration limit ({max_i} >= {_NITERMAX}) — "
        "divergence-free constraint not satisfied; check NITERMAX in BioReactor.c"
    )
