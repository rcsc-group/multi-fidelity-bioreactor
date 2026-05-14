"""Grid convergence: kLa_25 must agree within 20% between fidelity=5 and fidelity=6.

Physical basis: the first-order oxygen transfer coefficient kLa depends on
interface area and near-interface mass transfer.  At fidelity=5 (32×32 effective
cells) and fidelity=6 (64×64) the interface is already well-resolved for the
canonical geometry and forcing.  A difference > 20% signals that the fidelity=6
result has not yet been reached by fidelity=5 and the reference grid is under-
resolved.

Reference data (fidelity=6):
  runs/health_l6_video/ — t_end=100, canonical params, kLa_25=0.0832
  This directory is pre-computed and committed to runs/.  The test skips if it is
  absent so it does not break CI on a fresh clone without the reference data.

Fidelity=5 run (fidelity=5):
  Executed inline with t_end=55 (C* reaches 25% at t≈51.5 for canonical params).
  Marked hpc: intended to run on an OSCAR compute node, not on the login node.
  Typical runtime: ~20–40 min with OMP_NUM_THREADS=4.

Usage:
  pytest -m hpc tests/verification/test_grid_convergence.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[2]))
from tests.conftest import CANONICAL_PARAMS, run_bioreactor
from scripts.postprocess import main as postprocess_main

PROJECT_ROOT = Path(__file__).parents[2]
L6_REF_DIR   = PROJECT_ROOT / "runs" / "health_l6_video"
_KLA_RTOL    = 0.20   # 20 % relative tolerance between L5 and L6

# Run just past the 25 % saturation crossing (t≈51.5 for canonical params).
_PARAMS_L5 = {**CANONICAL_PARAMS, "run_id": "grid_conv_l5", "fidelity": 5, "t_end": 55.0}
_L5_TIMEOUT = 7200   # seconds — set for an OSCAR compute node with 4 threads


@pytest.mark.hpc
def test_kla_grid_converged_l5_vs_l6(tmp_path):
    """kLa_25 at fidelity=5 must be within 20 % of the fidelity=6 reference.

    Threshold calibration:
      L6 reference (health_l6_video): kLa_25 = 0.0832
      20 % tolerance leaves room for genuine grid-sensitivity while firmly
      catching a broken interface or Henry's-law term (which would give O(1)
      errors, not O(10 %) errors).
    """
    # ── load L6 reference ────────────────────────────────────────────────────
    results_json = L6_REF_DIR / "results.json"
    if not results_json.exists():
        pytest.skip(
            f"L6 reference data not found at {results_json}. "
            "Run the health_l6_video SLURM job first:\n"
            "  python scripts/simulate.py --fidelity 6 --t-end 100 --run-id health_l6_video"
        )
    l6 = json.loads(results_json.read_text())
    kla_l6 = l6.get("kLa_25")
    if kla_l6 is None or kla_l6 != kla_l6:   # None or NaN
        pytest.skip(f"L6 reference kLa_25 is not finite: {kla_l6}")

    # ── run fidelity=5 ───────────────────────────────────────────────────────
    run_dir = run_bioreactor(_PARAMS_L5, tmp_path, timeout=_L5_TIMEOUT)

    tr_oxy = run_dir / "tr_oxy.dat"
    if not tr_oxy.exists():
        pytest.fail(
            "tr_oxy.dat not written by fidelity=5 run — simulation may have crashed. "
            f"Check {run_dir} for stderr output."
        )

    results_l5 = postprocess_main(str(run_dir))
    kla_l5 = results_l5.get("kLa_25")

    if kla_l5 is None or kla_l5 != kla_l5:
        pytest.fail(
            f"kLa_25 is NaN for fidelity=5 run — C* may not have reached 25 % by "
            f"t_end={_PARAMS_L5['t_end']}.  Check tr_oxy.dat in {run_dir}."
        )

    # ── convergence assertion ────────────────────────────────────────────────
    rel_err = abs(kla_l5 - kla_l6) / kla_l6
    assert rel_err < _KLA_RTOL, (
        f"kLa_25 not grid-converged: L5={kla_l5:.4f}, L6={kla_l6:.4f}, "
        f"relative error={rel_err:.1%} (threshold {_KLA_RTOL:.0%}). "
        "Increase resolution or check interface/Henry's-law implementation."
    )
