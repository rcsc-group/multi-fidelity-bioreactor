"""Guard: a CONTINUATION restart must carry the tracer AND the oxygen across.

`restart_continue=1` (src/params_read.h) exists so one experiment can span a
walltime boundary: the soluble tracers keep their restored values and are not
re-injected. Figs 9/10 (mixing time) and 11/12 (kLa) all depend on it, and at
the time this file was written it had been verified exactly once, by hand, at
L5 -- with assertions on the TRACER only. Oxygen was observed in passing and
never asserted, and nothing in the suite guarded any of it.

The strongest available check is a split-vs-continuous comparison: run the
same physics once unbroken, and once as two segments joined by a continuation
restart, then require the two to agree at the end. That catches field loss,
re-injection, and anything that silently stops a source term at the seam --
without needing to know which of those went wrong.
"""
from __future__ import annotations

import math
import sys
import pathlib

import numpy as np
import pytest

PROJECT_ROOT = pathlib.Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from tests.conftest import CANONICAL_PARAMS, run_bioreactor, load_tr_oxy  # noqa: E402
from scripts.dump_fields import fields as dump_fields  # noqa: E402

# non-dimensional rocking period for CANONICAL_PARAMS (omega_b=3.93, 7 deg)
_L, _B = CANONICAL_PARAMS["geometry"]["a"], CANONICAL_PARAMS["geometry"]["b"]
_T_PER = 2 * math.pi / CANONICAL_PARAMS["omega_b"]
_H = 2 * _B
_V = _L / 4 * (_H + 0.5 * _L * math.tan(math.radians(CANONICAL_PARAMS["theta_max"][0])))
_U = _V / (_H * 0.5) / _T_PER
T_PER_ND = _T_PER / (_L / _U)

N_MIX = 2          # inject after 2 cycles
SEG = 3            # cycles per segment
_COL_OXY, _COL_C2 = 2, 8


def _params(run_id, t_end_cycles, **over):
    p = {**CANONICAL_PARAMS, "run_id": run_id, "fidelity": 4,
         "n_mix_cycles": N_MIX, "t_end": round(t_end_cycles * T_PER_ND, 6)}
    p.update(over)
    return p


def _tail(run_dir, col):
    a = load_tr_oxy(run_dir)
    assert a.shape[0] > 2, f"{run_dir.name}: tr_oxy.dat too short"
    return float(a[-1, col])


@pytest.mark.medium
def test_continuation_matches_an_unbroken_run(tmp_path):
    """Split (3+3 cycles, continuation) must equal continuous (6 cycles).

    Covers BOTH fields. The tracer is inert after injection, so agreement
    there tests only that the field survived the seam. Oxygen is REPLENISHED
    every timestep by `event oxygen (t=t_mix; i++)`, so agreement there also
    tests that the source term keeps firing after the restart -- on a
    continuation t_mix is recomputed relative to the new checkpoint and can
    land in the future, silently pausing replenishment.
    """
    cont = run_bioreactor(_params("cont", 2 * SEG), tmp_path, timeout=600)
    seg1 = run_bioreactor(_params("seg1", SEG), tmp_path, timeout=600)
    ck = seg1 / "checkpoint.dump"
    assert ck.exists(), "segment 1 wrote no checkpoint.dump"
    # t_checkpoint > 0 is what ARMS the restart (BioReactor.c gates the branch on
    # it, not on argv[2] being present). Without it the binary ignores the dump
    # and silently runs a fresh simulation -- which is how this test first
    # "failed", and how a production submit script was quietly cold-starting.
    t_ck = dump_fields(ck)[0]["t"]
    seg2 = run_bioreactor(_params("seg2", SEG, restart_continue=1, t_checkpoint=t_ck),
                          tmp_path, timeout=600, restart_from=ck)

    for name, col in (("tracer c2", _COL_C2), ("oxygen", _COL_OXY)):
        a, b = _tail(cont, col), _tail(seg2, col)
        assert a > 0, f"continuous run has no {name} at all -- test setup is wrong"
        assert b == pytest.approx(a, rel=0.10), (
            f"{name}: continuous={a:.6g} vs split-with-continuation={b:.6g} "
            f"(ratio {b/a:.3f}). A continuation must reproduce the unbroken "
            f"run. For oxygen, a LOW value means replenishment stopped at the "
            f"seam: `event oxygen (t=t_mix; i++)` with t_mix recomputed "
            f"relative to the new checkpoint does not fire until that future "
            f"time, so the gas side stops being topped up."
        )


@pytest.mark.medium
def test_default_restart_still_wipes_and_reinjects(tmp_path):
    """restart_continue=0 must keep chain.py's sweep semantics exactly.

    Each segment of a parameter sweep is a NEW tracer experiment, so the
    fields are zeroed and re-injected n_mix_cycles later. A fix aimed at
    continuations must not silently change this path.
    """
    seg1 = run_bioreactor(_params("d_seg1", SEG), tmp_path, timeout=600)
    ck = seg1 / "checkpoint.dump"
    assert ck.exists(), "segment 1 wrote no checkpoint.dump"
    t_ck = dump_fields(ck)[0]["t"]
    seg2 = run_bioreactor(_params("d_seg2", 1, t_checkpoint=t_ck),
                          tmp_path, timeout=600, restart_from=ck)

    a = load_tr_oxy(seg2)
    early = a[: max(2, len(a) // 4)]
    assert float(np.max(early[:, _COL_C2])) < 1e-9, (
        "restart_continue=0 must zero the tracer at the start of the segment "
        "(chain.py sweep semantics); found a non-zero value"
    )


@pytest.mark.medium
def test_kla_survives_a_restart_only_when_segments_are_stitched(tmp_path):
    """kLa across a continuation is correct -- but ONLY from the joined series.

    Measured: C* is continuous across the seam to 0.035%, and kLa fitted to the
    STITCHED series reproduces an unbroken run to ratio 1.000 at all three
    thresholds. So the solver side is sound.

    The trap is on the postprocessing side. postprocess.py computes kLa per RUN
    DIRECTORY, and a later segment restores an ALREADY-SATURATED oxygen field --
    so C* starts near 1 and crosses every threshold at its first row, returning
    the same kLa for 10%, 25% and 50%. That is the exact signature
    test_kla_values_differ_across_saturation_levels was written to catch, and a
    chained kLa sweep would hit it silently on every segment after the first.

    Guarded here because Figs 11/12 depend on it and chaining is how the long
    low-rpm points will have to run.
    """
    from scripts.postprocess import _compute_c_star, _kla_5pt_at_threshold

    NC = 20
    # require_complete on all three: this test compares runs against each
    # other, so a truncated one makes the comparison meaningless rather than
    # merely noisy. CI failed here from 2026-09-16 with continuous=1.057 vs
    # stitched=1.928 while the same test passed locally at ratio 1.000, and
    # the continuous run -- the only one of the three that is a single long
    # process -- is the one a slow runner would truncate first.
    cont = run_bioreactor(_params("k_cont", NC), tmp_path, timeout=1800,
                          require_complete=True)
    seg1 = run_bioreactor(_params("k_s1", NC // 2), tmp_path, timeout=1800,
                          require_complete=True)
    ck = seg1 / "checkpoint.dump"
    assert ck.exists(), "segment 1 wrote no checkpoint.dump"
    t_ck = dump_fields(ck)[0]["t"]
    seg2 = run_bioreactor(_params("k_s2", NC // 2, restart_continue=1,
                                  t_checkpoint=t_ck),
                          tmp_path, timeout=1800, restart_from=ck,
                          require_complete=True)

    tc, cc = _compute_c_star(cont)
    t1, c1 = _compute_c_star(seg1)
    t2, c2 = _compute_c_star(seg2)

    # (a) the oxygen field itself carries across
    assert abs(c2[0] - c1[-1]) < 0.01, (
        f"C* jumps {abs(c2[0]-c1[-1]):.4f} across the seam; the restored oxygen "
        f"field is not continuous")

    # (b) stitched kLa reproduces the unbroken run
    keep = t2 > t1[-1]
    ts = np.concatenate([t1, t2[keep]])
    cs = np.concatenate([c1, c2[keep]])
    order = np.argsort(ts)
    ts, cs = ts[order], cs[order]
    for thr in (0.10, 0.25, 0.50):
        a = _kla_5pt_at_threshold(tc, cc, thr)
        b = _kla_5pt_at_threshold(ts, cs, thr)
        if a != a:
            continue
        assert b == pytest.approx(a, rel=0.10), (
            f"kLa_{int(thr*100)}: continuous={a:.4g} vs stitched={b:.4g}\n"
            f"  continuous: {len(tc)} samples, t in [{tc[0]:.4f}, {tc[-1]:.4f}], "
            f"C* in [{cc.min():.4f}, {cc.max():.4f}]\n"
            f"  stitched:   {len(ts)} samples, t in [{ts[0]:.4f}, {ts[-1]:.4f}], "
            f"C* in [{cs.min():.4f}, {cs.max():.4f}]\n"
            f"  seam at t={t1[-1]:.4f}; if the two t-ranges differ the runs "
            f"are not the same experiment and the kLa gap is a symptom")

    # (c) the trap: a later segment ALONE must not be treated as a kLa run
    solo = [_kla_5pt_at_threshold(t2, c2, thr) for thr in (0.10, 0.25)]
    if all(v == v for v in solo):
        assert solo[0] == pytest.approx(solo[1], rel=1e-6), (
            "expected the known-degenerate case (a saturated segment crosses "
            "every threshold at its first row); if this no longer holds, "
            "per-segment kLa may have become meaningful and this guard needs "
            "rethinking rather than deleting")
