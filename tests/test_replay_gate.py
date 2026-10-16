"""Record/replay gate, pre-registered in diary.md 2026-10-16. The criteria below are fixed.

F0 no regression : the default binary before and after the change gives byte-identical tr_oxy.dat and
                   results.json on an L5 serial run.
F1 machinery     : L6, 32.5 rpm, theta 7. A coupled run records 3 periods from the release time; replaying
                   them from the start checkpoint gives oxy_liq_sum, c2_liq_sum and c2_liq_sum2 equal to
                   relative 1e-6 at every logged time (serial, one rank).
F2 looping       : L6, coupled run with release then 40 cycles; one recorded period replayed with looping for
                   the same 40 cycles. Reported, not asserted (the diary thresholds decide the physics):
                   max|C*_replay - C*_coupled|, kLa_1T_10 and dtmix_0.50 of both, and the period-to-period
                   flow defect. Asserted: the replay ran the whole 40 cycles and the numbers are finite.

These jobs run through scripts/simulate.submit_slurm (one rank, default account), so they are marked hpc
and slow. Run: uv run pytest tests/test_replay_gate.py -m "hpc and slow" -s
REPLAY_GATE_SRC points the F1/F2 builds at another source tree (the red step of TDD: the pre-change source).
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import replay_gate as rg  # noqa: E402

pytestmark = [pytest.mark.hpc, pytest.mark.slow]

BASE_REV = "b3a5f02"   # last commit before the record/replay change
SRC = Path(os.environ.get("REPLAY_GATE_SRC", ROOT / "src"))
F2_CYCLES = int(os.environ.get("REPLAY_GATE_F2_CYCLES", "40"))
F2_SPINUP = int(os.environ.get("REPLAY_GATE_F2_SPINUP", "40"))


def _base_src() -> Path:
    """Sources of BASE_REV, exported from git into the work area."""
    d = rg.WORK / f"src_{BASE_REV}"
    if not d.exists():
        d.mkdir(parents=True)
        names = subprocess.run(["git", "ls-tree", "--name-only", BASE_REV, "src/"], cwd=ROOT, capture_output=True,
                               text=True, check=True).stdout.split()
        for n in names:
            blob = subprocess.run(["git", "show", f"{BASE_REV}:{n}"], cwd=ROOT, capture_output=True, check=True)
            (d / Path(n).name).write_bytes(blob.stdout)
    return d


def _run(params: dict, checkpoint: Path | None = None, walltime: str = "01:00:00", src: Path = SRC) -> Path:
    job = rg.submit(params, checkpoint=checkpoint, walltime=walltime, src=src)
    state = rg.wait(job)
    assert state == "COMPLETED", f"job {job} ended {state}; see logs/slurm_{job}.err"
    return rg.run_dir(params["run_id"])


@pytest.fixture(scope="module")
def binaries():
    """default (current src), record and replay builds of the source under test; all fresh."""
    out = {}
    for tag, flags in (("default", []), ("rec", ["-DRECORD_FLUX"]), ("replay", ["-DREPLAY_FLUX"])):
        out[tag] = rg.build_binary(f"{tag}-test", flags, src=SRC)
        assert rg.binary_is_current(out[tag], SRC), f"{out[tag]} is older than {SRC}"
    return out


def test_f0_default_build_is_unchanged(binaries):
    base = rg.build_binary(f"base-{BASE_REV}", [], src=_base_src())
    runs = {}
    for tag, binary in (("before", base), ("after", binaries["default"])):
        p = rg.make_params(f"f0_{tag}", binary, fidelity=5, n_mix_cycles=2, post_cycles=3)
        runs[tag] = _run(p, walltime="00:30:00", src=_base_src() if tag == "before" else SRC)
    for name in ("tr_oxy.dat", "results.json", "vol_frac_interf.dat", "normf.dat", "shear_stress.dat"):
        a, b = (runs["before"] / name).read_bytes(), (runs["after"] / name).read_bytes()
        assert a == b, f"{name} differs between the pre-change and the post-change default binary"


def _f1_runs(binaries):
    rec_p = rg.make_params("f1_coupled", binaries["rec"], fidelity=6, n_mix_cycles=10, post_cycles=3,
                           n_rec_periods=3, rec_dir=str(rg.scratch_dir("f1_coupled")))
    coupled = _run(rec_p)
    t_per = rg.t_period(rec_p)
    for name in ("flux_rec.0.bin", "replay_start.dump"):
        assert (rg.scratch_dir("f1_coupled") / name).exists(), f"the coupled run wrote no {name}"
    rep_p = rg.make_params("f1_replay", binaries["replay"], fidelity=6, n_mix_cycles=0, post_cycles=0,
                           t_checkpoint=10 * t_per, restart_continue=1, n_rec_periods=3, replay_loop=0,
                           rec_dir=str(rg.scratch_dir("f1_coupled")))
    rep_p["t_end"] = 3 * t_per
    replay = _run(rep_p, checkpoint=rg.scratch_dir("f1_coupled") / "replay_start.dump")
    return coupled, replay, t_per


def test_f1_machinery_exact(binaries):
    coupled, replay, t_per = _f1_runs(binaries)
    a, b = rg.load_tr_oxy(coupled / "tr_oxy.dat"), rg.load_tr_oxy(replay / "tr_oxy.dat")
    cols = [2, 8, 9]   # oxy_liq_sum, c2_liq_sum, c2_liq_sum2
    d = rg.rel_diff_by_time(a, b, cols)
    print("\nF1: common times", d["n_common"], "replay rows", d["n_b"], "max rel diff", d["max"])
    t_rel = d["t_common"]
    assert d["n_common"] == d["n_b"] > 3 * t_per / 0.02 * 0.9, "replay rows missing from the coupled log"
    # the replayed window must be the 3 periods after the release, not a trivial all-zero comparison
    assert np.max(np.abs(b[:, 8])) > 0.0 and np.max(np.abs(b[:, 2])) > 0.0
    for c in cols:
        assert d["max"][c] <= 1e-6, f"column {c}: max relative difference {d['max'][c]:.3e} > 1e-6"
    # the last row is the last output time before the window end (no row at the end itself: the run stops first)
    assert 13 * t_per - 0.02 - 1e-6 <= t_rel[-1] < 13 * t_per


def test_f2_looping_reports(binaries):
    rec_dir = rg.scratch_dir("f2_coupled")
    rec_p = rg.make_params("f2_coupled", binaries["rec"], fidelity=6, n_mix_cycles=F2_SPINUP,
                           post_cycles=F2_CYCLES, n_rec_periods=F2_CYCLES, rec_dir=str(rec_dir))
    t_per = rg.t_period(rec_p)
    reuse = os.environ.get("REPLAY_GATE_REUSE") == "1" and (rg.run_dir("f2_coupled") / "results.json").exists()
    coupled = rg.run_dir("f2_coupled") if reuse else _run(rec_p, walltime="04:00:00")
    rep_p = rg.make_params("f2_replay", binaries["replay"], fidelity=6, n_mix_cycles=0, post_cycles=0,
                           t_checkpoint=F2_SPINUP * t_per, restart_continue=1, n_rec_periods=1, replay_loop=1, rec_dir=str(rec_dir))
    rep_p["t_end"] = F2_CYCLES * t_per
    replay = _run(rep_p, checkpoint=rec_dir / "replay_start.dump", walltime="02:00:00")

    tc, cc = rg.cstar_series(coupled)
    tr, cr = rg.cstar_series(replay)
    m = tr >= F2_SPINUP * t_per
    dc = np.abs(cr[m] - np.interp(tr[m], tc, cc))
    res_c, res_r = rg.load_results(coupled), rg.load_results(replay)
    defect = rg.flow_defect(rec_dir, t_per, F2_CYCLES)
    out = {"max_dCstar": float(dc.max()), "t_at_max": float(tr[m][dc.argmax()]),
           "kLa_1T_10": {"coupled": res_c["kLa_1T_10"], "replay": res_r["kLa_1T_10"]},
           "dtmix_0.50": {"coupled": res_c["dtmix_0.50"], "replay": res_r["dtmix_0.50"]},
           "flow_defect_max": defect["max"], "flow_defect_per_period": defect["per_period"], "umax": defect["umax"],
           "replay_t_end": float(tr[-1]), "n_periods": F2_CYCLES}
    (rg.WORK / "f2_metrics.json").write_text(json.dumps(out, indent=2))
    print("\nF2:", json.dumps(out, indent=2))
    assert tr[-1] >= (F2_SPINUP + F2_CYCLES) * t_per * (1 - 1e-6) - t_per * 0.02, "replay stopped early"
    assert all(math.isfinite(x) for x in (out["max_dCstar"], out["flow_defect_max"]))
