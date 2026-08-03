"""Physical verification: the simulated fluid domain area matches geometry.a/b.

geometry.b is documented (docs_site/reference/params.md) as a HALF-height
semi-axis -- the embedded solid is built as |y| < b_nd, so the full bag
height is 2*b_nd and the total fluid-domain area (within the L0=1 box) is
L0 * 2*b_nd = 2*(geometry.b/geometry.a).

This is a direct structural check on the actual constructed geometry (the
solid()/y_fill construction, which was NOT the source of the 2026-08-03
bug -- that bug was in H_bio's formula and in geometry.b's default VALUE,
both orthogonal to whether solid() faithfully builds |y|<b_nd from
whatever b_nd currently is). This test is a general regression guard for
the geometry construction itself, complementing (not duplicating)
test_kim_fig_a16_velocity_rms.py, which is the test that actually would
have caught the 2026-08-03 bug (it checks against Kim et al.'s published
values, not against the code's own internal formulas).

Failure modes caught:
  - Any future off-by-factor-of-N error introduced while refactoring the
    solid()/y_fill construction (e.g. reverting toward upstream's
    0.5*Ly-style convention without updating b_nd's own definition).
  - A geometry.a/b unit or ordering mix-up.
"""
import sys
import pathlib

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))
from tests.conftest import CANONICAL_PARAMS, load_normf, run_bioreactor


@pytest.mark.medium
def test_fluid_domain_area_matches_geometry_a_b(tmp_path):
    """normf.dat's *_liq_vol columns (total fluid-domain area) must equal
    L0 * 2*(geometry.b/geometry.a) to within 2% -- L0=1 always (BioReactor.c
    non-dimensionalizes length by L_bio=geometry.a, so the domain box is
    always exactly 1 unit wide in these units).
    """
    params = {**CANONICAL_PARAMS, "run_id": "bag_height_check", "t_end": 1.0}
    run_dir = run_bioreactor(params, tmp_path)

    data = load_normf(run_dir)
    assert len(data) >= 2, "Too few output rows — simulation may not have run"

    # Omega_liq_vol (col 4), ux_liq_vol (col 8), uy_liq_vol (col 12) are all
    # the same normf() "volume" reduction (total fluid-domain area) —
    # check the last row (VOF hasn't moved much by t_end=1.0, but the
    # embedded-solid area itself is constant throughout the run anyway).
    measured_area = data[-1, 4]

    geometry_a = params["geometry"]["a"]
    geometry_b = params["geometry"]["b"]
    expected_area = 1.0 * 2 * (geometry_b / geometry_a)  # L0 * 2*b_nd

    rel_err = abs(measured_area - expected_area) / expected_area
    assert rel_err < 0.02, (
        f"Simulated fluid-domain area is {measured_area:.4f}, expected "
        f"{expected_area:.4f} (2*geometry.b/geometry.a) — off by "
        f"{rel_err:.1%} (measured/expected={measured_area/expected_area:.3f}). "
        "The embedded solid()/y_fill construction no longer matches "
        "geometry.b's documented half-height meaning."
    )
