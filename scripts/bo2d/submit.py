"""Submissions for the DFD BO comparison (2026-11-07): L9-only BO vs L7-L8-L9 chained BO.

cold(level, rpm, th): cold start, the L8 campaign protocol (kmix_l8_rpm32.5 params: 80-cycle spin-up, release,
150 cycles) at any level. Records every job in scripts/bo2d/state/jobs.jsonl with the arm it belongs to.
Real sbatch only with go=True (tests block sbatch through tests/conftest.py).
"""
from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))
import campaign as C  # noqa: E402

STATE = ROOT / "scripts" / "bo2d" / "state"
JOBS = STATE / "jobs.jsonl"
RANKS = {7: 8, 8: 16, 9: 32}
N_SPIN = 80
WALL = {7: "06:00:00", 8: "16:00:00", 9: "16:00:00"}
TAG = "bo"


def cold_params(level: int, rpm: float, th: float, tag: str = TAG, spin_only: bool = False) -> dict:
    p = json.loads((ROOT / "runs/kmix_l8_rpm32.5/params.json").read_text())
    for k in ("_canonical_run_dir", "next_run_id", "_err_file", "_source_run"):
        p.pop(k, None)
    p.update(run_id=C.run_id(level, rpm, th, tag=tag), fidelity=int(level), omega_b=2 * math.pi * rpm / 60.0,
             theta_max=[float(th), 0.0, 0.0], _binary=C.BINARY_L8)
    if spin_only:   # start state for a warm start only: N_SPIN cycles, ends at the release (no scalar data used)
        from scripts.postprocess import _t_scales
        _, T = _t_scales(p)
        p.update(t_end=round(N_SPIN * T, 6), n_mix_cycles=N_SPIN)
    return p


def record(**kw) -> None:
    STATE.mkdir(parents=True, exist_ok=True)
    with JOBS.open("a") as f:
        f.write(json.dumps(dict(kw, submitted=datetime.now(timezone.utc).isoformat())) + "\n")


def cold(level: int, rpm: float, th: float, arm: str, go: bool = False, tag: str = TAG,
         spin_only: bool = False) -> str | None:
    from scripts import simulate
    p = cold_params(level, rpm, th, tag=tag, spin_only=spin_only)
    if (ROOT / "runs" / p["run_id"]).exists():
        print(f"skip {p['run_id']}: run dir exists"); return None
    if not go:
        print(f"dry: cold L{level} {p['run_id']} ranks {RANKS[level]} arm {arm}"); return None
    jid = simulate.submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=WALL[level],
                                ntasks=RANKS[level], exclude=C.EXCLUDE, partition="batch", account="default",
                                cpus=1, mem="4G")
    record(job=str(jid), level=level, run=p["run_id"], rpm=rpm, theta=th, ranks=RANKS[level], arm=arm,
           start="spin" if spin_only else "cold")
    return str(jid)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("level", type=int); ap.add_argument("rpm", type=float); ap.add_argument("theta", type=float)
    ap.add_argument("--arm", default="seed"); ap.add_argument("--go", action="store_true")
    a = ap.parse_args()
    print(cold(a.level, a.rpm, a.theta, a.arm, a.go))


# ---- warm starts (2026-11-07). Same forcing as the source only: a theta change at restart is NOT physical in this
# solver (BioReactor.c l.186-189: w_bio_st = w_bio*T_bio depends on theta_max, so sin(w_bio_st*t) jumps in phase at
# t_checkpoint ~ 140; the u rescale at l.360 uses omega_b_prev/omega_b only). So: cross-LEVEL warm starts at the same
# (rpm, theta) only.
N_WARM, N_POST = 30, 40      # 30 settle cycles before release (diary 2026-11-07: L9 wave ~27 cycles)
BIN_XLEVEL = C.BINARY_L9


def warm_params(level: int, source: str, rpm: float, th: float, tag: str = TAG) -> dict:
    from scripts.postprocess import _t_scales
    sp = json.loads((ROOT / "runs" / source / "params.json").read_text())
    if abs(float(sp["omega_b"]) * 60 / (2 * math.pi) - rpm) > 0.01 or abs(float(sp["theta_max"][0]) - th) > 0.01:
        raise ValueError(f"{source}: forcing differs from ({rpm}, {th}); cross-condition restart is not physical")
    p = dict(sp)
    for k in ("_canonical_run_dir", "next_run_id", "_err_file", "_source_run"):
        p.pop(k, None)
    _, T = _t_scales(p)
    p.update(run_id=C.run_id(level, rpm, th, tag=tag), fidelity=int(level), t_checkpoint=float(C.source_checkpoint_time(source)),
             n_mix_cycles=N_WARM, t_end=round((N_WARM + N_POST) * T, 6), _binary=BIN_XLEVEL, _source_run=source)
    return p


def warm(level: int, source: str, rpm: float, th: float, arm: str, go: bool = False, tag: str = TAG):
    from scripts import simulate
    p = warm_params(level, source, rpm, th, tag=tag)
    if (ROOT / "runs" / p["run_id"]).exists():
        print(f"skip {p['run_id']}: run dir exists"); return None
    sd = ROOT / "runs" / source
    if not go:
        print(f"dry: warm L{level} {p['run_id']} from {source} arm {arm}"); return None
    jid = simulate.submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=WALL[level],
                                ntasks=RANKS[level], template=ROOT / "config" / "slurm_mpi_template.sh",
                                checkpoint=str((sd / "checkpoint.dump").resolve()), exclude=C.EXCLUDE,
                                partition="batch", account="default", cpus=1, mem="4G")
    rp = ROOT / "runs" / p["run_id"] / "params.json"
    rp.write_text(json.dumps(C.warm_provenance(json.loads(rp.read_text()), source, ROOT / "logs" / f"slurm_{jid}.err"), indent=1))
    record(job=str(jid), level=level, run=p["run_id"], rpm=rpm, theta=th, ranks=RANKS[level], arm=arm, start="warm",
           source=source)
    return str(jid)
