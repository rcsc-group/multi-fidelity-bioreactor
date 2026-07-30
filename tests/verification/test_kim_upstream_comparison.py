"""Informational-only comparison: our fork vs a near-literal reproduction of
Kim et al.'s own upstream driver code, at matched resolution and condition.

This test NEVER asserts and NEVER fails on a numeric mismatch. It exists
purely to surface, via a pytest warning, how far this project's fork has
drifted from Kim et al.'s own published code -- so a large, unexpected
*change* in that drift (e.g. introduced by a future refactor) is visible
to a human running it, instead of going unnoticed.

CI status: marked `medium`, DOES run in GitHub Actions' `medium-tests` job
-- Basilisk is already built from source there by `.github/actions/setup-
basilisk` (no OSCAR dependency; `_find_qcc()` below checks `$BASILISK`,
which CI's setup-basilisk step sets, before falling back to the OSCAR
persistent-storage install for local/manual runs).

Why this is warning-only, not a pass/fail gate
-----------------------------------------------
As of 2026-07-30 (diary.md), our fork's serial baseline and a minimal-diff
reproduction of Kim's own code (`tests/fixtures/kim_upstream/`, see its
README for exact provenance) disagree by ~1.74x in ux_rms/U_b at matched
condition and resolution (fidelity 6 / NN=64) -- for real, intentional,
already-documented reasons: a much shorter forcing ramp (3 rocking cycles
vs. Kim's own fixed 30-second ramp), a superellipse tank cross-section vs.
a literal rectangle, and a 2x difference in the liquid-volume convention
used to normalize normf() output. None of that is a bug in the ordinary
sense -- it's this project's fork intentionally diverging from Kim's setup
for its own optimization pipeline. A hard-tolerance assert here would
either be so loose it catches nothing, or so tight it fails today, on
main, for reasons unrelated to code health. Use test_mpi_checkpoint_parity.py
for the actual pass/fail regression guards.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[2]))
from tests.conftest import CANONICAL_PARAMS, load_normf, run_bioreactor

PROJECT_ROOT = Path(__file__).parents[2]
FIXTURE_DIR  = Path(__file__).parent.parent / "fixtures" / "kim_upstream"

# Same condition as CANONICAL_PARAMS, at fidelity 6 (NN=64) -- matches this
# fixture's hardcoded `const double NN = 64` in BioReactor.c exactly, so the
# two runs are resolution-matched without needing to patch the vendored file.
# t_end=15 gives ~8 rocking cycles of post-ramp data even under Kim's own
# (much longer) 30-second ramp -- see _kim_t_ramp_nd() below.
_OUR_PARAMS = {**CANONICAL_PARAMS, "run_id": "kim_cmp_our_fork", "fidelity": 6, "t_end": 15.0}
_L_BIO, _ANGLE_DEG, _RPM = 0.25, 7.0, 32.5   # argv for the vendored kim_upstream binary


def _find_qcc() -> Path | None:
    """Same discovery order as the project Makefile: $BASILISK env var (set by
    CI's setup-basilisk action, or by a developer's shell) first, then the
    OSCAR persistent-storage install, then whatever's on PATH.
    """
    import os
    import shutil as _shutil
    basilisk_env = os.environ.get("BASILISK")
    if basilisk_env and (Path(basilisk_env) / "qcc").exists():
        return Path(basilisk_env) / "qcc"
    oscar_qcc = Path("/oscar/data/dharri15/eaguerov/basilisk/src/qcc")
    if oscar_qcc.exists():
        return oscar_qcc
    found = _shutil.which("qcc")
    return Path(found) if found else None


@pytest.fixture(scope="module")
def kim_upstream_binary(tmp_path_factory):
    """Compile the vendored Kim upstream fixture with whichever qcc is available."""
    qcc = _find_qcc()
    if qcc is None:
        pytest.skip("no qcc found ($BASILISK unset, no OSCAR install, none on PATH)")
    # qcc gets confused mixing an absolute source path with a different cwd
    # (matches a known qcc quirk: it wants sources resolvable relative to its
    # own working directory) -- copy the fixture into the build dir first and
    # compile with a plain relative filename instead.
    import shutil as _shutil
    build_dir = tmp_path_factory.mktemp("kim_upstream_build")
    for f in FIXTURE_DIR.glob("*.[ch]"):
        _shutil.copy(f, build_dir / f.name)
    basilisk_src = qcc.parent
    binary = build_dir / "BioReactor_kim_upstream"
    result = subprocess.run(
        [str(qcc), "-O2", "-Wall", "-disable-dimensions",
         "BioReactor.c", "-o", str(binary),
         "-L", str(basilisk_src / "gl"), "-lglutils", "-lfb_tiny",
         "-lm", "-lGLU", "-lGL", "-lX11"],
        cwd=build_dir, capture_output=True, text=True,
        env={**__import__("os").environ, "BASILISK": str(basilisk_src)},
    )
    if result.returncode != 0 or not binary.exists():
        pytest.skip(f"Kim upstream fixture failed to compile:\n{result.stderr[-2000:]}")
    return binary


def _run_kim_upstream(binary: Path, run_dir: Path, timeout: int) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    for sub in ("Data_all", "Data_specific", "Fig_oxy", "Fig_tr", "Fig_vol", "Fig_vor"):
        (run_dir / sub).mkdir(exist_ok=True)
    try:
        subprocess.run(
            [str(binary.resolve()), str(_L_BIO), str(_ANGLE_DEG), str(_RPM)],
            cwd=run_dir, capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        pass


def _kim_t_ramp_nd() -> float:
    """Kim's own fixed 30-second ramp, in the nondimensional time units used
    by normf.dat -- NOT the fork's 3-rocking-cycle ramp (see test_grid_
    convergence._t_ramp_nd, which is fork-specific and would be wrong here).
    """
    import math as m
    omega_b = 2 * m.pi * _RPM / 60
    T_per   = 2 * m.pi / omega_b
    th      = m.radians(_ANGLE_DEG)
    H       = 0.0715  # matches this fixture's `Ly=0.286` * L_bio=0.25
    V       = _L_BIO / 4 * (H + 0.5 * _L_BIO * m.tan(th))
    U_bio   = V / (H * 0.5) / T_per
    T_bio   = _L_BIO / U_bio
    return 30.0 / T_bio


def _mean_post_ramp_vel_rms(normf_data: np.ndarray, t_ramp: float) -> float:
    t       = normf_data[:, 1]
    vel_rms = np.sqrt(normf_data[:, 7] ** 2 + normf_data[:, 11] ** 2)
    post    = vel_rms[t > t_ramp]
    if len(post) < 5:
        return float("nan")
    return float(post[len(post) // 2:].mean())


@pytest.mark.medium
def test_our_fork_vs_kim_upstream_informational(kim_upstream_binary, tmp_path):
    """Report (never assert) the velocity-RMS gap between our fork and Kim's
    own upstream code, at matched condition and resolution."""
    our_dir = run_bioreactor(_OUR_PARAMS, tmp_path, timeout=1800)
    if not (our_dir / "normf.dat").exists():
        pytest.skip("our fork's normf.dat not written -- run failed, nothing to compare")

    kim_dir = tmp_path / "kim_cmp_upstream"
    _run_kim_upstream(kim_upstream_binary, kim_dir, timeout=1800)
    if not (kim_dir / "normf.dat").exists():
        pytest.skip("Kim upstream fixture's normf.dat not written -- run failed")

    from tests.verification.test_grid_convergence import _t_ramp_nd as _our_t_ramp_nd
    vel_our = _mean_post_ramp_vel_rms(load_normf(our_dir), _our_t_ramp_nd(_OUR_PARAMS))
    vel_kim = _mean_post_ramp_vel_rms(load_normf(kim_dir), _kim_t_ramp_nd())

    if math.isnan(vel_our) or math.isnan(vel_kim):
        pytest.skip(f"insufficient post-ramp data (our={vel_our}, kim={vel_kim})")

    ratio = vel_our / vel_kim if vel_kim else float("nan")
    warnings.warn(
        f"our fork vs Kim upstream (informational, not a failure): "
        f"post-ramp combined velocity RMS -- our fork={vel_our:.4f}, "
        f"Kim upstream={vel_kim:.4f}, ratio={ratio:.3f}. "
        "Historical reference (2026-07-30, diary.md): ratio ~0.57 (5.42/9.44 "
        "in peak ux_rms/U_b terms, a related but not identical metric). "
        "A large *change* from that historical reference is worth investigating; "
        "the discrepancy itself is expected and not a bug -- see this file's docstring."
    )
