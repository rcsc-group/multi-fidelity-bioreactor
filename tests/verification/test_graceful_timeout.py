"""Guard: a walltime kill must never destroy a run's work.

`event dump_checkpoint (t = t_dump_checkpoint)` fires exactly ONCE, at t_end.
A run killed by SLURM before reaching it leaves nothing. On 2026-09-21
fig9_l9_rpm20 was 49% through 24 h of L9 compute with no checkpoint on disk,
and fig9_l10_seg1 18 h in with only the seed dump it had restarted FROM; both
were on track to lose everything, and neither could be rescued because SLURM
refuses to raise TimeLimit on a RUNNING job.

MECHANISM. Basilisk ships this in src/maxruntime.h: an `i += 10` event that
all-reduces perf.t with MPI_MAX and dumps once the wall clock comes within
five minutes of a declared limit. A deadline beats catching SIGUSR1 -- no
async-signal-safety constraints, no dependence on SLURM delivering a signal
to every MPI rank at a comparable timestep, and MPI_MAX gives rank consensus
for free, where a signal needs an explicit agreement step or one rank enters
collective dump() alone and hangs the job it was trying to save.

THE HARD PART IS PHASE, NOT THE DUMP. Every checkpoint written so far lands
on a zero-crossing of its writer's forcing, and the restart path exploits
that -- it INFERS the writer's phase by snapping to the nearest period
boundary:

    t_phase_offset = T_per_st*round (t_ck/T_per_st) - t_ck   (BioReactor.c)

A deadline dump arrives at whatever phase the clock ran out at. Snapping it
moves the forcing by up to half a period, so the resumed segment starts with
the bag at a different angle than the one it stopped at -- the "escape"
failure documented at BioReactor.c:470 and measured on a theta 2 -> 7 -> 7
chain. So the writer must RECORD its phase instead of letting the reader
guess. Both halves are pinned below.
"""
from __future__ import annotations

import json
import math
import pathlib
import re
import sys
import time

import numpy as np
import pytest

PROJECT_ROOT = pathlib.Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from tests.conftest import CANONICAL_PARAMS, run_bioreactor  # noqa: E402

SRC = PROJECT_ROOT / "src" / "BioReactor.c"
TEMPLATE = PROJECT_ROOT / "config" / "slurm_mpi_template.sh"

_L, _B = CANONICAL_PARAMS["geometry"]["a"], CANONICAL_PARAMS["geometry"]["b"]
_T_PER = 2 * math.pi / CANONICAL_PARAMS["omega_b"]
_H = 2 * _B
_V = _L / 4 * (_H + 0.5 * _L * math.tan(
    math.radians(CANONICAL_PARAMS["theta_max"][0])))
_U = _V / (_H * 0.5) / _T_PER
T_PER_ND = _T_PER / (_L / _U)

PHASE_SIDECAR = "checkpoint.phase"
ENV_DEADLINE = "BIOREACTOR_MAXRUNTIME_S"
# Must match MAXRUNTIME_RESERVE_S in src/BioReactor.c.
RESERVE_S = 300
# Must match N_RAMP_CYCLES in src/BioReactor.c.
N_RAMP_CYCLES = 3


# --------------------------------------------------------------- source-level
# Millisecond tests. Most of the risk is a missing piece rather than a wrong
# one: until every part is wired the feature fails silently as "no dump, as
# before", which looks identical to the bug it fixes.

def test_solver_has_a_walltime_deadline_event():
    src = SRC.read_text()
    assert ENV_DEADLINE in src, "solver never learns its own walltime"


def test_deadline_is_agreed_across_ranks():
    """dump() is collective: every rank must decide to stop on the same step.

    Without the reduction one rank can cross the deadline a timestep before
    the others, enter dump() alone, and hang the job until the very kill the
    feature exists to pre-empt -- strictly worse than not having it.
    """
    src = SRC.read_text()
    assert re.search(r"mpi_all_reduce\s*\(\s*[\w.]+\s*,\s*MPI_DOUBLE\s*,"
                     r"\s*MPI_MAX\s*\)", src), \
        "no MPI_MAX reduction guarding the deadline check"


def test_emergency_dump_records_its_phase():
    """The writer must write its phase down, not leave it to be inferred."""
    src = SRC.read_text()
    assert PHASE_SIDECAR in src, "no phase sidecar written"


def test_restart_prefers_a_recorded_phase_over_guessing():
    """The reader must use the sidecar when present.

    The nearest-boundary form stays as the fallback so every checkpoint
    written before this feature still restores correctly.
    """
    src = SRC.read_text()
    assert "round (t_ck/T_per_st)" in src, \
        "the legacy inference was removed; old checkpoints would break"
    i_side = src.index(PHASE_SIDECAR)
    assert src.count(PHASE_SIDECAR) >= 2, \
        "sidecar is written but never read back"
    assert i_side > 0


def test_slurm_template_tells_the_job_its_deadline():
    txt = TEMPLATE.read_text()
    assert ENV_DEADLINE in txt, "template never exports the deadline"


def test_reserve_exceeds_the_time_a_dump_takes():
    """The reserve must cover writing a checkpoint at the largest level.

    An L10 checkpoint is ~196 MB (measured on fig9_l10_seg1). A reserve
    shorter than the write leaves a truncated dump, which is worse than none
    because it looks resumable.
    """
    src = SRC.read_text()
    m = re.search(r"MAXRUNTIME_RESERVE_S\s+(\d+)", src)
    assert m, "no named reserve constant"
    assert int(m.group(1)) >= 300, "reserve too short for a 196 MB write"


# --------------------------------------------------------------- behavioural

@pytest.mark.medium
def test_deadline_produces_a_resumable_checkpoint(tmp_path):
    """Interrupt a run by deadline; get a checkpoint back.

    t_end is far beyond what the deadline allows, so the run CANNOT reach
    dump_checkpoint on its own: any checkpoint on disk is the deadline's
    doing.
    """
    params = {**CANONICAL_PARAMS, "run_id": "deadline", "fidelity": 4,
              "n_mix_cycles": 1, "t_end": round(200 * T_PER_ND, 6)}
    t0 = time.time()
    run_dir = run_bioreactor(params, tmp_path, timeout=600,
                             env={ENV_DEADLINE: "330"})
    elapsed = time.time() - t0

    ck = run_dir / "checkpoint.dump"
    assert ck.exists() and ck.stat().st_size > 0, \
        "deadline did not produce a checkpoint"
    assert (run_dir / PHASE_SIDECAR).exists(), \
        "no phase sidecar: the restart would have to guess the phase"
    assert elapsed < 300, \
        f"ran {elapsed:.0f}s against a 330s deadline with a 300s reserve"


@pytest.mark.medium
def test_emergency_checkpoint_preserves_phase(tmp_path):
    """Forcing must be continuous across an emergency restart.

    Run the same physics unbroken, and interrupted mid-period then resumed.
    The bag angle is a pure function of t and the forcing phase, so if the
    resumed segment reconstructs the phase the two theta(t) agree; if it
    snaps to the nearest boundary they differ by up to half a period, which
    at these amplitudes is a sign flip.

    Compared on the SIGNED spatial-mean u_x (normf.dat, ux_liq_savg) rather
    than on a mixing scalar. The forcing drives the sloshing directly, so
    u_x oscillates at the rocking period and in phase with it; a half-period
    error flips its sign, which is the largest signal any logged quantity
    offers. A rectified quantity (rms, or normf's |.|-averaged columns) would
    be blind to exactly that error, and chi would show it only indirectly and
    many cycles later.
    """
    # Long enough that a usable window remains AFTER the ramp: the cut lands
    # near cycle 4, the ramp runs 3 more, and everything past cycle ~8 is
    # fully driven. Measured 100% sign agreement over 185 samples at these
    # settings; at 8 cycles the whole comparison sat inside the ramp and read
    # 71%, which was the observable failing, not the feature.
    cycles = 16.0
    base = {**CANONICAL_PARAMS, "fidelity": 4, "n_mix_cycles": 1}

    whole = {**base, "run_id": "phase_whole",
             "t_end": round(cycles * T_PER_ND, 6)}
    t0 = time.time()
    whole_dir = run_bioreactor(whole, tmp_path, timeout=900)
    whole_s = time.time() - t0
    ref = _forcing_series(whole_dir)

    # Cut the second run partway through by giving it a deadline it hits
    # mid-run. Derived from the reference's measured cost rather than
    # hard-coded: the same fraction holds whatever the machine's speed, and a
    # fixed number would either let the run finish (no interruption, nothing
    # tested) or stop it before the ramp ends.
    deadline = RESERVE_S + max(3.0, 0.30 * whole_s)
    part = {**base, "run_id": "phase_part",
            "t_end": round(cycles * T_PER_ND, 6)}
    part_dir = run_bioreactor(part, tmp_path, timeout=900,
                              env={ENV_DEADLINE: f"{deadline:.0f}"})
    assert (part_dir / "checkpoint.dump").exists(), "no emergency checkpoint"
    assert (part_dir / PHASE_SIDECAR).exists(), "no phase sidecar"
    side = json.loads((part_dir / PHASE_SIDECAR).read_text())
    t_ck = float(side["t"])
    # The whole point: this dump is NOT on a period boundary.
    off = abs((t_ck / T_PER_ND) % 1.0 - 0.5)
    assert off < 0.49, ("emergency dump landed on a period boundary by luck; "
                        "the adversarial case was not exercised")

    resumed = {**base, "run_id": "phase_resume",
               "t_checkpoint": t_ck, "restart_continue": 1,
               "t_end": round((cycles * T_PER_ND) - t_ck, 6)}
    res_dir = run_bioreactor(
        resumed, tmp_path, timeout=900,
        restart_from=part_dir / "checkpoint.dump",
        extra_files=[part_dir / PHASE_SIDECAR])

    got = _forcing_series(res_dir)
    ts = got[:, 0]
    # Compared on SIGN, and only AFTER the ramp.
    #
    # A continuation ramps the forcing amplitude 0 -> 1 over
    # N_RAMP_CYCLES = 3 cycles (BioReactor.c:439). Sign is amplitude-blind,
    # which is why it beats comparing values -- but u_x is the FLUID's
    # response, not the forcing, and under a ramped-down forcing the
    # response's phase lag shifts too. So inside the ramp even the sign
    # disagrees for reasons that have nothing to do with the checkpoint
    # (measured: 71-77% agreement, entirely within the ramp). Past it the
    # forcing is at full amplitude and the comparison means what it claims.
    swing = float(np.ptp(ref[:, 1]))
    assert swing > 0, "reference forcing signal is flat; nothing to compare"
    want = np.interp(ts, ref[:, 0], ref[:, 1])
    # Past the ramp, and away from crossings where either sign is noise.
    m = (ts > t_ck + (N_RAMP_CYCLES + 0.5) * T_PER_ND) & \
        (np.abs(want) > 0.2 * swing)
    assert m.sum() > 50, f"only {m.sum()} post-ramp samples to compare"
    agree = float(np.mean(np.sign(got[m, 1]) == np.sign(want[m])))

    # Negative control. 100% agreement is only evidence if this metric can
    # fail: shift the reference by half a period -- precisely what a
    # nearest-boundary phase snap does -- and require the same comparison to
    # collapse. Without this, an insensitive metric would read as a pass.
    shifted = np.interp(ts + 0.5 * T_PER_ND, ref[:, 0], ref[:, 1])
    ctrl = float(np.mean(np.sign(got[m, 1]) == np.sign(shifted[m])))
    assert ctrl < 0.35, (
        f"negative control failed: a deliberate half-period shift still "
        f"agrees {ctrl*100:.0f}% of the time, so this comparison cannot "
        f"detect the error it exists to detect")

    assert agree > 0.9, (
        f"forcing phase not preserved across the emergency restart: sign "
        f"agrees on {agree*100:.0f}% of {m.sum()} post-ramp samples "
        f"(half-period control: {ctrl*100:.0f}%).")


# --------------------------------------------------------------- helpers

# normf.dat column order (BioReactor.c header):
#   0 i  1 t  2..5 Omega_*  6..9 ux_*  10..13 uy_*  14 ux_liq_savg  ...
_COL_T, _COL_UX_SAVG = 1, 14


def _forcing_series(run_dir: pathlib.Path) -> np.ndarray:
    """(t, ux_liq_savg) -- the signed mean velocity, which tracks the forcing.

    Signed is essential. normf's ux_liq_avg is normf(v).avg, which sums |v|
    and so is identically blind to the sign flip a half-period phase error
    produces; ux_liq_savg is the dv()-weighted signed mean added for Fig 2
    precisely because the rectified column hid this class of error.
    """
    f = run_dir / "normf.dat"
    if not f.exists():
        pytest.skip(f"{run_dir.name}: no normf.dat")
    a = np.loadtxt(f, skiprows=1)
    if a.ndim != 2 or a.shape[1] <= _COL_UX_SAVG:
        pytest.skip(f"{run_dir.name}: normf.dat has no signed columns "
                    f"(binary predates 2026-09-16)")
    return a[:, [_COL_T, _COL_UX_SAVG]]
