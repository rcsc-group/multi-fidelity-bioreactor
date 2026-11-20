"""Feature tests at L9 on the HEAD solver (diary 2026-11-20, pre-registered X9 and R9).

X9  cross-level warm start: L8 cold source (80 cycles) -> L9 warm (12 warm-up, release, 40 post; campaign protocol)
    vs L9 cold from rest (80-cycle spin-up, release, 40 post). Sites (25, 4) and (32.5, 7).
R9  replay: from the settled L9 warm run, a coupled RECORD run (release 1 cycle after the restart, the first K = 2
    periods recorded, 40 cycles coupled) and a REPLAY run looping those K periods for 40 cycles. Sites (32.5, 7) and
    (37.5, 7); (37.5, 7) gets its own L8 source and L9 warm run.

Stages (run one at a time; each stage needs the previous one's runs finished):
    build    compile the HEAD binaries (MPI lean + cross-level; MPI lean + record; MPI lean + replay) on a compute node
    sources  L8 sources (3 sites) and L9 cold runs (2 sites)
    warm     L9 warm runs from the finished L8 sources (3 sites)
    record   R9 record runs from the finished L9 warm runs
    replay   R9 replay runs from the finished record runs
    dry      print what each stage would submit
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402
from scripts.campaign2d import campaign as cp  # noqa: E402
from scripts.simulate import submit_slurm  # noqa: E402

HEAD = "ee811a9"   # src commit the binaries were built from (src unchanged since; later commits touch tools/docs only)
import os
SMOKE = os.environ.get("FT9_SMOKE") == "1"   # the same chain at L5 -> L6, short, 4 ranks: checks the driver only
TAG = "ft9s" if SMOKE else "ft9"
BIN = {k: rg.BIN_DIR / f"BioReactor-mpi-{k}-{HEAD}" for k in ("xlevel", "record", "replay")}
FLAGS = {"xlevel": ["-DEXTRA_TRACERS=0", "-DCROSS_LEVEL_WARMSTART=1"],
         "record": ["-DEXTRA_TRACERS=0", "-DCROSS_LEVEL_WARMSTART=1", "-DRECORD_FLUX=1"],
         "replay": ["-DEXTRA_TRACERS=0", "-DREPLAY_FLUX=1"]}
X9 = [(25.0, 4.0), (32.5, 7.0)]
R9 = [(32.5, 7.0), (37.5, 7.0)]
SITES = sorted(set(X9 + R9))
N_L8, N_SPIN, N_WARM, N_POST, K = (4, 4, 12, 3, 2) if SMOKE else (80, 80, 12, 40, 2)
LO, HI = (5, 6) if SMOKE else (8, 9)
EXCLUDE = "node2356,node2410,node2411"
RUNS = ROOT / "runs"
TEMPLATE = ROOT / "config" / "slurm_mpi_template.sh"


def rid(kind: str, rpm: float, th: float) -> str:
    return f"{TAG}_{kind}_rpm{rpm:g}_th{th:g}"


def base(level: int, rpm: float, th: float, binary: Path, run_id: str) -> dict:
    """Production parameters (a c2d_l8 run) at this level, site and binary."""
    p = json.loads((RUNS / "c2d_l8_rpm32.19_th4.51" / "params.json").read_text())
    p = {k: v for k, v in p.items() if not k.startswith("_") and not k.endswith("_prev")
         and k not in ("t_checkpoint", "restart_continue", "next_run_id")}
    p.update(run_id=run_id, fidelity=level, omega_b=2 * 3.141592653589793 * rpm / 60.0,
             theta_max=[th, 0.0, 0.0], _binary=str(binary))
    return p


def T(p: dict) -> float:
    return pp._t_scales(p)[1]


def submit(p: dict, ranks: int, wall: str, checkpoint: Path | None = None, dry: bool = False) -> str:
    if dry:
        print("  ", p["run_id"], p["fidelity"], f"t_end={p['t_end']:.3f}", "ck" if checkpoint else "cold", Path(p["_binary"]).name)
        return "dry"
    if (RUNS / p["run_id"]).exists():
        raise SystemExit(f"{p['run_id']} exists: refusing to overwrite")
    if SMOKE:
        ranks, wall = 4, "01:00:00"
    j = submit_slurm(p, project_root=ROOT, runs_root=RUNS, walltime=wall, ntasks=ranks, template=TEMPLATE, cpus=1,
                     mem="4G", checkpoint=str(checkpoint) if checkpoint else None, exclude=EXCLUDE, partition="batch")
    print("  ", p["run_id"], "job", j)
    return j


def last_time(run: str) -> float:
    return cp.source_checkpoint_time(run)


def stage_build(dry: bool) -> None:
    for k, b in BIN.items():
        print(k, b, *FLAGS[k])
        if not dry:
            rg.build_binary(b.name.removeprefix("BioReactor-"), FLAGS[k], mpi=True)


def stage_sources(dry: bool) -> None:
    for b in BIN.values():
        if not dry and not b.exists():
            raise SystemExit(f"missing {b}: run the build stage")
    for rpm, th in SITES:
        p = base(LO, rpm, th, BIN["xlevel"], rid("l8", rpm, th))
        p.update(n_mix_cycles=N_L8, t_end=round(N_L8 * T(p), 6))
        submit(p, 16, "16:00:00", dry=dry)
    for rpm, th in X9:
        p = base(HI, rpm, th, BIN["xlevel"], rid("l9cold", rpm, th))
        p.update(n_mix_cycles=N_SPIN, t_end=round((N_SPIN + N_POST) * T(p), 6))
        submit(p, 32, "48:00:00", dry=dry)


def stage_warm(dry: bool) -> None:
    for rpm, th in SITES:
        src = rid("l8", rpm, th)
        sp = json.loads((RUNS / src / "params.json").read_text())
        p = cp.l9_warm_params(sp, src, 0.0 if dry else last_time(src), rpm, th, str(BIN["xlevel"]))
        p.update(run_id=rid("l9warm", rpm, th), fidelity=HI, n_mix_cycles=N_WARM if not SMOKE else 1,
                 t_end=round(((N_WARM if not SMOKE else 1) + N_POST) * T(p), 6))
        submit(p, 32, "24:00:00", RUNS / src / "checkpoint.dump", dry=dry)


def stage_record(dry: bool) -> None:
    for rpm, th in R9:
        src = rid("l9warm", rpm, th)
        sp = json.loads((RUNS / src / "params.json").read_text())
        p = {k: v for k, v in sp.items() if not k.startswith("_")}
        p.update(run_id=rid("rec", rpm, th), t_checkpoint=0.0 if dry else last_time(src), n_mix_cycles=1,
                 t_end=round((1 + N_POST) * T(sp), 6), n_rec_periods=K, rec_dir=".", _binary=str(BIN["record"]))
        p.pop("restart_continue", None)
        submit(p, 32, "24:00:00", RUNS / src / "checkpoint.dump", dry=dry)


def stage_replay(dry: bool) -> None:
    for rpm, th in R9:
        src = rid("rec", rpm, th)
        sp = json.loads((RUNS / src / "params.json").read_text())
        rdir = RUNS / src
        p = {k: v for k, v in sp.items() if not k.startswith("_") and k not in ("t_end", "n_mix_cycles")}
        p.update(run_id=rid("rep", rpm, th), n_mix_cycles=0, t_end=round(N_POST * T(sp), 6), restart_continue=1,
                 n_rec_periods=K, replay_loop=1, rec_dir=str(rdir), replay_scalar_level=HI, replay_diag=0,
                 replay_project_flux=0, _binary=str(BIN["replay"]))
        submit(p, 32, "12:00:00", rdir / "replay_start.dump", dry=dry)


STAGES = {"build": stage_build, "sources": stage_sources, "warm": stage_warm, "record": stage_record,
          "replay": stage_replay}

if __name__ == "__main__":
    st = sys.argv[1]
    if st == "dry":
        for name, f in STAGES.items():
            print("==", name)
            try:
                f(True)
            except FileNotFoundError as e:
                print("   (needs earlier stage)", e.filename)
    else:
        STAGES[st](False)
