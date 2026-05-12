"""Verify that geometry.a is wired to L_bio (the characteristic length).

L_bio is used for all non-dimensionalization. If it is hardcoded at 0.25 m,
changing geometry.a has no effect on the dynamics: two runs with different
geometry.a produce identical non-dim output. If L_bio is correctly wired,
the non-dim time scale T_bio = L_bio/U_bio differs → posY_max oscillates at
a different non-dim rate between the two runs.
"""
import sys
import pathlib
import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))
from tests.conftest import CANONICAL_PARAMS, run_bioreactor, load_vol_frac


@pytest.mark.medium
def test_geometry_a_changes_nondim_dynamics(tmp_path):
    """posY_max must differ between geometry.a=0.20 and geometry.a=0.25 runs.

    If L_bio is hardcoded, both runs use the same non-dimensionalization and
    produce identical dynamics. If wired, T_bio ∝ f(L_bio) differs and the
    interface traces diverge.
    """
    params_020 = {
        **CANONICAL_PARAMS,
        "run_id": "geo_a_020",
        "geometry": {"a": 0.20, "b": 0.071, "n": 8.0},
        "t_end": 3.0,
    }
    params_025 = {
        **CANONICAL_PARAMS,
        "run_id": "geo_a_025",
        "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
        "t_end": 3.0,
    }

    run_dir_020 = run_bioreactor(params_020, tmp_path)
    run_dir_025 = run_bioreactor(params_025, tmp_path)

    vf_020 = load_vol_frac(run_dir_020)
    vf_025 = load_vol_frac(run_dir_025)

    posY_020 = vf_020[:, 4]  # posY_max column
    posY_025 = vf_025[:, 4]

    n = min(len(posY_020), len(posY_025))
    assert n >= 4, "Too few output rows to compare — check t_end or output frequency"

    mean_diff = np.abs(posY_020[:n] - posY_025[:n]).mean()
    assert mean_diff > 1e-4, (
        f"posY_max is identical for geometry.a=0.20 vs 0.25 (mean diff={mean_diff:.2e}): "
        "L_bio is likely hardcoded and not wired to geometry.a"
    )
