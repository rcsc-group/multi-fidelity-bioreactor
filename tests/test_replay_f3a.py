"""F3a (diary.md 2026-10-16, PRE-REGISTERED test F3): scalars and f on a finer grid than the flow.

L6 flow, 5 periods recorded after the release (same case as F1), replayed at scalar levels Ls = 6, 7, 8 with the
recorded uf prolonged by refine_face_solenoidal. Measured per Ls and written to WORK/f3a_ls<Ls>.json:
  div_far_over_umax / div_all_over_umax : max |div uf| Delta / max|uf| (cells with cs = 1 and four open faces / all)
  liquid_volume_drift, c2_mass_drift    : max over the replay of |x(t) - x(0)| / x(0), from vol_frac_interf.dat / tr_oxy.dat
Ls = 6 must equal the coupled run (as F1, relative 1e-6). The numeric thresholds of the diary (div <= 1e-10, drifts
<= 1e-6) are reported, and asserted at the end by test_f3a_summary, so a failure is a result, not a test bug.
Run: uv run pytest tests/test_replay_f3a.py -m "hpc and slow" -s
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import replay_gate as rg  # noqa: E402

pytestmark = [pytest.mark.hpc, pytest.mark.slow]
PERIODS, SPINUP, CID = 5, 10, "f3a_coupled"


def _run(params, checkpoint=None, walltime="04:00:00"):
    job = rg.submit(params, checkpoint=checkpoint, walltime=walltime)
    state = rg.wait(job)
    assert state == "COMPLETED", f"job {job} ended {state}; see logs/slurm_{job}.err"
    return rg.run_dir(params["run_id"])


@pytest.fixture(scope="module")
def binaries():
    return {"rec": rg.build_binary("rec-f3", ["-DRECORD_FLUX"]), "replay": rg.build_binary("replay-f3", ["-DREPLAY_FLUX"])}


@pytest.fixture(scope="module")
def coupled(binaries):
    p = rg.make_params(CID, binaries["rec"], fidelity=6, n_mix_cycles=SPINUP, post_cycles=PERIODS,
                       n_rec_periods=PERIODS, rec_dir=str(rg.scratch_dir(CID)))
    return _run(p)


@pytest.mark.parametrize("ls", [6, 7, 8])
def test_f3a_scalar_mesh(binaries, coupled, ls):
    cp = json.loads((coupled / "params.json").read_text())
    t_per = rg.t_period(cp)
    p = rg.make_params(f"f3a_replay_ls{ls}", binaries["replay"], fidelity=6, n_mix_cycles=0, post_cycles=0,
                       t_checkpoint=SPINUP * t_per, restart_continue=1, n_rec_periods=PERIODS, replay_loop=0,
                       rec_dir=str(rg.scratch_dir(CID)), replay_scalar_level=ls, replay_diag=1)
    p["t_end"] = PERIODS * t_per
    replay = _run(p, checkpoint=rg.scratch_dir(CID) / "replay_start.dump")
    diag = np.loadtxt(rg.scratch_dir(f"f3a_replay_ls{ls}") / "replay_diag.txt")
    far, allc, umax = diag[:, 2], diag[:, 3], diag[:, 4]
    tr, vf = rg.load_tr_oxy(replay / "tr_oxy.dat"), rg.load_tr_oxy(replay / "vol_frac_interf.dat")
    out = {"ls": ls, "div_far_over_umax": float(far.max() / umax.max()), "div_all_over_umax": float(allc.max() / umax.max()),
           "umax": float(umax.max()), "liquid_volume_drift": float(np.max(np.abs(vf[:, 2] - vf[0, 2])) / vf[0, 2]),
           "c2_mass_drift": float(np.max(np.abs(tr[:, 8] - tr[0, 8])) / abs(tr[0, 8])), "n_rows": len(tr)}
    if ls == 6:
        d = rg.rel_diff_by_time(rg.load_tr_oxy(coupled / "tr_oxy.dat"), tr, [2, 8, 9])
        out["vs_coupled_max_rel"], out["n_common"] = d["max"], d["n_common"]
    (rg.WORK / f"f3a_ls{ls}.json").write_text(json.dumps(out, indent=2))
    print("\nF3a:", json.dumps(out))
    if ls == 6:
        assert d["n_common"] == len(tr) and all(v <= 1e-6 for v in d["max"].values())


def test_f3a_summary():
    """The diary thresholds, applied to what the three replays measured."""
    res = [json.loads((rg.WORK / f"f3a_ls{ls}.json").read_text()) for ls in (6, 7, 8)]
    for r in res[1:]:
        assert r["div_far_over_umax"] <= 1e-10, r
    for r in res:
        assert r["liquid_volume_drift"] <= 1e-6 and r["c2_mass_drift"] <= 1e-6, r
