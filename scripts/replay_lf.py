"""Low-fidelity oxygen-transfer estimate by looping ONE recorded flow period (diary 2026-10-21, RESULT periodic-flow replay).

Per condition (rpm, theta) at L8, two chained Slurm jobs (16 ranks each, default account, node2356 excluded):
  RECORD  restart from the source run's checkpoint.dump (same level and condition; *_prev = current forcing),
          restart_continue = 0 and n_mix_cycles = 1: the restart zeroes tracer and oxygen, releases them 1 cycle
          after the checkpoint (t_mix = t_ck + 1 T), and the RECORD_FLUX binary records exactly one period from the
          release (t_rec0 unset -> t_mix; n_rec_periods = 1). Run length 2 T.
  REPLAY  BioReactor-replayL16 restores the record's replay_start.dump (restart_continue = 1, n_mix_cycles = 0),
          loops the one recorded period (replay_loop = 1, n_rec_periods = 1) for 22 cycles on the L8 scalar grid.
The replay job starts only when the record job ended OK (--dependency=afterok).
Source runs: kmix_l8_rpm<r> (theta 7) and t6_l8_rpm<r>_th<theta> (theta 4, 5.5, 9).
"""
from __future__ import annotations

import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import replay_gate as rg  # noqa: E402
from scripts import postprocess as pp  # noqa: E402
from scripts.simulate import submit_slurm  # noqa: E402

SRC_RUNS = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs")
REC_ROOT = Path("/oscar/data/dharri15/eaguerov/replay_records/lf")
REC_BIN = rg.BIN_DIR / "BioReactor-recL"
REP_BIN = rg.BIN_DIR / "BioReactor-replayL16"
SRC_TREE = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/src")  # sources the binaries were built from
RPMS = [15 + 2.5 * i for i in range(10)]
THETAS = [7.0, 4.0, 5.5, 9.0]
N_REPLAY_CYCLES = 22
N_MIX = 1
NTASKS = 16
EXCLUDE = "node2356"
SOURCE_CYCLE = 230.0


def source_id(rpm: float, theta: float) -> str:
    return f"kmix_l8_rpm{rpm:g}" if theta == 7.0 else f"t6_l8_rpm{rpm:g}_th{theta:g}"


def lf_id(rpm: float, theta: float) -> str:
    return f"lf_l8_rpm{rpm:g}_th{theta:g}"


def rec_id(rpm: float, theta: float) -> str:
    return f"lfrec_l8_rpm{rpm:g}_th{theta:g}"


def checkpoint_time(dump: Path) -> float:
    """Absolute non-dimensional time stored in the first 8 bytes of a Basilisk dump (what the solver reads)."""
    return float(np.fromfile(dump, dtype="<f8", count=1)[0])


def record_params(src_params: dict, rid: str, t_ck: float, rec_dir: Path) -> dict:
    zeros = [0.0, 0.0, 0.0]
    _, T = pp._t_scales(src_params)
    p = dict(src_params)
    p.update(run_id=rid, t_checkpoint=t_ck, t_end=round(2 * T, 6), n_mix_cycles=N_MIX,
             omega_b_prev=src_params["omega_b"], theta_max_prev=list(src_params["theta_max"]), phi_angular_prev=zeros,
             amplitude_h_prev=zeros, phi_horizontal_prev=zeros, omega_h_prev=0.0, _binary=str(REC_BIN),
             n_rec_periods=1, rec_dir=str(rec_dir))
    p.pop("restart_continue", None)   # = 0: tracer and oxygen zeroed, released at t_ck + N_MIX T
    p.pop("t_rec0", None)             # unset: the record starts at the release time
    return p


def replay_params(src_params: dict, rid: str, t_rec0: float, rec_dir: Path) -> dict:
    """Same keys as scripts/replay_gate_f3/drv_f3b_loop.py (t_checkpoint is overwritten by the dump's own time)."""
    _, T = pp._t_scales(src_params)
    p = {k: v for k, v in src_params.items() if k not in ("t_end", "n_mix_cycles")}
    p.update(run_id=rid, n_mix_cycles=0, t_checkpoint=t_rec0, t_end=N_REPLAY_CYCLES * T, restart_continue=1,
             n_rec_periods=1, replay_loop=1, rec_dir=str(rec_dir), replay_scalar_level=8, replay_diag=0,
             replay_project_flux=0, _binary=str(REP_BIN))
    return p


def conditions() -> list[tuple[float, float]]:
    return [(r, th) for th in THETAS for r in RPMS]


def _submit(params: dict, checkpoint: Path, walltime: str, dependency: str | None = None) -> str:
    for d in (rg.run_dir(params["run_id"]), rg.scratch_dir(params["run_id"])):
        if d.exists():
            shutil.rmtree(d)
    rg.RUNS.mkdir(parents=True, exist_ok=True)
    return submit_slurm(params, project_root=ROOT, runs_root=rg.RUNS, walltime=walltime,
                        template=ROOT / "config" / "slurm_mpi_template.sh", cpus=1, ntasks=NTASKS, mem="4G",
                        checkpoint=str(checkpoint), dependency=dependency, exclude=EXCLUDE,
                        mpi_scratch_root=rg.MPI_SCRATCH)


def submit_pair(rpm: float, theta: float, record_only: bool = False, dry_run: bool = False) -> dict:
    """Submit the record job and (afterok) the replay job of one condition. Returns the ids."""
    for b in (REC_BIN, REP_BIN):
        if not rg.binary_is_current(b, SRC_TREE):
            raise RuntimeError(f"{b} is older than {SRC_TREE}: rebuild before submitting")
    src = SRC_RUNS / source_id(rpm, theta)
    ck = src / "checkpoint.dump"
    sp = json.loads((src / "params.json").read_text())
    _, T = pp._t_scales(sp)
    t_ck = checkpoint_time(ck)
    cyc = t_ck / T
    if abs(cyc - SOURCE_CYCLE) > 0.01:
        raise RuntimeError(f"{src.name}: checkpoint at cycle {cyc:.4f}, expected {SOURCE_CYCLE:g}")
    rid, lid = rec_id(rpm, theta), lf_id(rpm, theta)
    rdir = REC_ROOT / lid
    info = dict(rpm=rpm, theta=theta, source=src.name, cycle=cyc, rec_run=rid, rep_run=lid, rec_dir=str(rdir))
    if dry_run:
        return info
    if rdir.exists():
        shutil.rmtree(rdir)
    rdir.mkdir(parents=True)
    pr = record_params(sp, rid, t_ck, rdir)
    job_rec = _submit(pr, ck, "02:00:00")
    info["job_record"] = job_rec
    if record_only:
        return info
    t_rec0 = t_ck + N_MIX * T
    pp_ = replay_params(sp, lid, t_rec0, rdir)
    # submit_slurm copies the checkpoint when it submits, and the record's replay_start.dump does not exist yet: stage
    # the source dump as a placeholder, then point the staged file at the record's dump. The replay cannot start
    # before the record ended OK (afterok), so the link is in place in time.
    job_rep = _submit(pp_, ck, "02:00:00", dependency=f"afterok:{job_rec}")
    staged = rg.scratch_dir(lid) / "checkpoint.dump"
    staged.unlink()
    staged.symlink_to(rg.scratch_dir(rid) / "replay_start.dump")
    info["job_replay"] = job_rep
    return info


QOIS = ("kLa_1T_10", "kLa_1T_25", "kLa_1T_50", "dtmix_0.50", "dtmix_0.75")


def control_params(src_params: dict, rid: str, t_ck: float, cycles: float = 23.0) -> dict:
    """Coupled continuation with the record's protocol (restart at the checkpoint, release 1 cycle later), source binary,
    no recording: the reference that separates 'the replay is wrong' from 'the flow at cycle 231 differs from cycle 81'."""
    p = record_params(src_params, rid, t_ck, Path("."))
    _, T = pp._t_scales(src_params)
    for k in ("n_rec_periods", "rec_dir"):
        p.pop(k)
    p.update(t_end=round(cycles * T, 6), _binary=src_params["_binary"])
    return p
