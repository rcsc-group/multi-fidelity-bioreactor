"""Regression guards for the CROSS-LEVEL WARM START (converge at a coarse
fidelity, refine() onto a finer grid, continue there).

Why these tests exist
----------------------
2026-09-14 (diary.md): the cross-level warm start was silently producing
wrong physics for four days. `fs` (the embedded-boundary face fractions) does
not survive dump/restore -- Basilisk's `dump_list()` skips face fields
outright (src/output.h:1037, `if (!s.face && ...)`) while `cm`/`cs` is a
plain scalar and round-trips perfectly. What came back was not zero but
STRUCTURALLY broken: roughly half of all FULLY OPEN faces (both neighbouring
cells cs > 0.99, no geometry anywhere near) were hard ZERO.

`refine_embed_linear` (src/embed-tree.h:268) gates every interpolation branch
on fs truthiness -- and the "pathological cases" fallback's 1D-gradient
corrections are fs-gated too, leaving `s[] = coarse(s)`. So during refine()
about half the cells fell to piecewise-constant injection while the
interleaved rest interpolated properly: a period-2 checkerboard imposed over
the WHOLE domain. Cost: +29.7% (L6->L7) and +42.2% (L7->L8) on domain-mean
wall shear, and +219.8% / +126.4% on viscous dissipation -- with no crash, no
warning, and a plausible-looking flow field. Eight fixes aimed at the
near-wall region failed because the metric used to judge them (tau) is itself
near-wall-only; dissipation, a VOLUME quantity, stayed pegged at +198% and
flat for 58 cycles and was what finally exposed the true extent.

Fix: call solid() on the coarse grid BEFORE refine(). The two tests here
guard the two independent things that went wrong:

  test_fs_is_sane_after_restart   -- the INVARIANT. A face whose both
      neighbours are entirely fluid cannot have zero open area. This is the
      cheap check that would have caught the bug on day one; the C code now
      asserts it on every restart (assert_fs_sane in src/BioReactor.c) and
      aborts, so this test asserts the run does NOT abort.

  test_cross_level_matches_cold_start -- the OUTCOME. A warm start refined
      from a coarser checkpoint must converge to the same state as a cold
      start at the fine fidelity, on BOTH a wall quantity (tau) and a volume
      quantity (dissipation). Dissipation is the one that matters: it is
      ~|grad u|^2 and is therefore the sensitive detector of grid-scale
      noise that wall metrics miss. Never drop it from this comparison.

CI status: marked `medium`, serial, fidelity 4->5, no MPI or SLURM needed.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[2]))
from tests.conftest import CANONICAL_PARAMS, load_normf, run_bioreactor

PROJECT_ROOT = Path(__file__).parents[2]

_COARSE_FIDELITY = 4
_FINE_FIDELITY = 5
_SPINUP_T = 8.0    # coarse run: long enough to get past the 3-cycle ramp
_COMPARE_T = 6.0   # forward run at fine fidelity, warm and cold alike
_TIMEOUT = 1800

# Measured on the real L6->L7 and L7->L8 pilots (diary.md 2026-09-14): with
# the fix both pairs land at +0.01-0.04% on tau and +0.01-0.17% on
# dissipation, versus +29.7%/+219.8% and +42.2%/+126.4% broken. At the tiny
# fidelities used here the two runs are not expected to agree that tightly,
# but the broken state is off by TENS of percent, so a loose threshold still
# separates them decisively without being flaky.
_RTOL_TAU = 0.15
_RTOL_EDISS = 0.30


@pytest.fixture(scope="module")
def crosslevel_binary():
    """Build BioReactor-crosslevel (-DCROSS_LEVEL_WARMSTART=1) via the Makefile."""
    build = subprocess.run(
        ["make", "build-crosslevel"], cwd=PROJECT_ROOT, capture_output=True, text=True,
    )
    binary = PROJECT_ROOT / "build" / "BioReactor-crosslevel"
    if build.returncode != 0 or not binary.exists():
        pytest.skip(f"build-crosslevel failed:\n{build.stderr}")
    return binary


def _coarse_checkpoint(tmp_path):
    """Run a short coarse simulation and return (its dir, its final time)."""
    params = {
        **CANONICAL_PARAMS,
        "run_id": "xlevel_coarse_src",
        "fidelity": _COARSE_FIDELITY,
        "t_end": _SPINUP_T,
    }
    run_dir = run_bioreactor(params, tmp_path, timeout=_TIMEOUT)
    dump = run_dir / "checkpoint.dump"
    if not dump.exists():
        pytest.fail(f"coarse run wrote no checkpoint.dump in {run_dir}")
    return run_dir, float(load_normf(run_dir)[-1, 1])


def _run_warmstart(binary, tmp_path, dump, t_checkpoint, run_id):
    """Restart `dump` (coarse) at the FINE fidelity -> triggers refine()."""
    params = {
        **CANONICAL_PARAMS,
        "run_id": run_id,
        "fidelity": _FINE_FIDELITY,
        "t_end": _COMPARE_T,
        "t_checkpoint": t_checkpoint,
        "omega_b_prev": CANONICAL_PARAMS["omega_b"],
        "theta_max_prev": CANONICAL_PARAMS["theta_max"],
    }
    run_dir = tmp_path / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "params.json").write_text(json.dumps(params))
    shutil.copy(dump, run_dir / "checkpoint.dump")
    proc = subprocess.run(
        [str(binary.resolve()), "params.json", "checkpoint.dump"],
        cwd=run_dir, capture_output=True, text=True, timeout=_TIMEOUT,
    )
    return run_dir, proc


def _tau_and_ediss(run_dir: Path):
    """Time-averaged domain-mean wall shear and dissipation over the last half."""
    dat = np.loadtxt(run_dir / "shear_stress.dat", skiprows=1)
    if dat.ndim != 2 or len(dat) < 8:
        pytest.fail(f"shear_stress.dat in {run_dir} too short -- run likely crashed")
    half = len(dat) // 2
    return float(np.sqrt(np.mean(dat[half:, 5] ** 2))), float(np.sqrt(np.mean(dat[half:, 9] ** 2)))


@pytest.mark.medium
def test_fs_is_sane_after_restart(crosslevel_binary, tmp_path):
    """A fully-open face (both neighbours cs>0.99) can never have fs==0.

    src/BioReactor.c asserts this on every restart and exits(1) with a
    diagnostic if violated, so a clean exit IS the assertion. Guards against
    anyone reordering solid() back after refine().
    """
    src_dir, t_ck = _coarse_checkpoint(tmp_path)
    run_dir, proc = _run_warmstart(
        crosslevel_binary, tmp_path, src_dir / "checkpoint.dump", t_ck, "xlevel_fs_sane"
    )
    assert "fs is not usable" not in (proc.stderr or ""), (
        "fs failed its sanity invariant after restart -- face fractions are "
        "stale/partial where the geometry is nowhere near. This is the "
        "2026-09-14 bug: solid() must recompute cs/fs BEFORE any refine() "
        f"that interpolates u/p/g.\n\nstderr:\n{proc.stderr[-2000:]}"
    )
    assert proc.returncode == 0, (
        f"cross-level warm start exited {proc.returncode}\n"
        f"stderr:\n{(proc.stderr or '')[-2000:]}"
    )


@pytest.mark.medium
def test_cross_level_matches_cold_start(crosslevel_binary, tmp_path):
    """A cross-level warm start must reach the same state as a cold start.

    Checks a wall quantity (tau) AND a volume quantity (dissipation). The
    volume one is not optional: the 2026-09-14 bug left tau ~30% off but
    dissipation ~200% off, and near-wall-only metrics are structurally blind
    to a domain-wide defect.
    """
    src_dir, t_ck = _coarse_checkpoint(tmp_path)
    warm_dir, proc = _run_warmstart(
        crosslevel_binary, tmp_path, src_dir / "checkpoint.dump", t_ck, "xlevel_warm"
    )
    assert proc.returncode == 0, (
        f"warm start exited {proc.returncode}\nstderr:\n{(proc.stderr or '')[-2000:]}"
    )

    cold_params = {
        **CANONICAL_PARAMS,
        "run_id": "xlevel_cold",
        "fidelity": _FINE_FIDELITY,
        "t_end": t_ck + _COMPARE_T,
    }
    cold_dir = run_bioreactor(cold_params, tmp_path, timeout=_TIMEOUT)

    tau_w, ed_w = _tau_and_ediss(warm_dir)
    tau_c, ed_c = _tau_and_ediss(cold_dir)
    tau_err = abs(tau_w - tau_c) / (tau_c + 1e-30)
    ed_err = abs(ed_w - ed_c) / (ed_c + 1e-30)

    assert tau_err < _RTOL_TAU, (
        f"cross-level warm start disagrees with cold start on wall shear: "
        f"warm={tau_w:.6g}, cold={tau_c:.6g}, rel err={tau_err:.1%} "
        f"(threshold {_RTOL_TAU:.0%}). See diary.md 2026-09-14."
    )
    assert ed_err < _RTOL_EDISS, (
        f"cross-level warm start disagrees with cold start on DISSIPATION: "
        f"warm={ed_w:.6g}, cold={ed_c:.6g}, rel err={ed_err:.1%} "
        f"(threshold {_RTOL_EDISS:.0%}). Dissipation ~|grad u|^2 is the "
        "sensitive detector of grid-scale noise that wall metrics miss -- in "
        "the 2026-09-14 bug it was +198% while tau showed only +30%. "
        "Do not relax this threshold without reading that diary entry."
    )
