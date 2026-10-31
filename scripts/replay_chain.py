"""Chained mfbml with the offline replay in place of the L8/L9 DNS (diary 2026-10-31).
Per condition of the database run (theta 7: kmix_l8_rpm*, fig9_l9_rpm*): a recording DNS from rest to the release
(cycle 80) plus ONE recorded cycle, then a replay of that cycle looped for the rest of the database window.
Same binaries as the pre-registered loop test (pr4-rec / pr4-rep: old seed, old sums, as the database runs).
  diff LEVEL RPM    params of our recording run against the database run (must differ only in the replay keys)
  dns LEVEL [RPM..] recording DNS
  replay LEVEL      afterok of the DNS: the looped replay to the database t_end
  report LEVEL      kLa (C* = 0.10) and dtmix 0.50 of replay against database, per rpm
Run ids: rc{LEVEL}_rpm{RPM}_{rec|rep}.
Usage: uv run python scripts/replay_chain.py diff 8 25 | dns 8 | replay 8 | report 8
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402
from scripts import prereg_offline as po  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402
from scripts.test_h1_mfbml import RPMS  # noqa: E402

THETA = 7.0
DB = {8: "kmix_l8_rpm{:g}", 9: "fig9_l9_rpm{:g}"}
REC_ROOT = Path("/oscar/scratch/eaguerov/replay_records/chain")   # ~4 GB/cycle at L8, ~17 at L9
REPLAY_KEYS = {"run_id", "_binary", "n_rec_periods", "rec_dir", "t_end", "n_mix_cycles", "t_checkpoint",
               "restart_continue", "replay_scalar_level", "replay_diag", "replay_project_flux", "replay_loop"}


def rid(level: int, rpm: float, kind: str) -> str:
    return f"rc{level}_rpm{rpm:g}_{kind}"


def db_params(level: int, rpm: float) -> dict:
    return json.loads((ROOT / "runs" / DB[level].format(rpm) / "params.json").read_text())


def rec_params(level: int, rpm: float) -> dict:
    """The database params, with the recording binary, stopped one cycle after the release."""
    p = dict(db_params(level, rpm))
    _, T = pp._t_scales(p)
    d = REC_ROOT / rid(level, rpm, "rec")
    p.update(run_id=rid(level, rpm, "rec"), _binary=str(po.BIN["rec"]), n_rec_periods=1, rec_dir=str(d),
             t_end=(p["n_mix_cycles"] + 1.05) * T)
    return p


def diff(level: int, rpm: float) -> None:
    a = db_params(level, rpm)
    b = rg.make_params("x", po.BIN["rec"], fidelity=level, n_mix_cycles=a["n_mix_cycles"], post_cycles=1.0,
                       rpm=rpm, theta=THETA)
    for k in sorted(set(a) | set(b)):
        if k not in REPLAY_KEYS and a.get(k) != b.get(k):
            print(f"{k}: database {a.get(k)!r}  make_params {b.get(k)!r}")


def dns(level: int, rpms: list[float]) -> None:
    f = rg.WORK / f"replay_chain_dns_{level}.json"
    jobs = json.loads(f.read_text()) if f.exists() else {}
    for rpm in rpms:
        p = rec_params(level, rpm)
        d = Path(p["rec_dir"])
        d.mkdir(parents=True, exist_ok=True)
        jobs[f"{rpm:g}"] = po._submit(p, po.NT_A, "24:00:00" if level >= 9 else "06:00:00")
    f.write_text(json.dumps(jobs))
    print(jobs)


def replay(level: int) -> None:
    jobs = json.loads((rg.WORK / f"replay_chain_dns_{level}.json").read_text())
    out = {}
    for key, jid in jobs.items():
        rpm = float(key)
        src = rec_params(level, rpm)
        _, T = pp._t_scales(src)
        db = db_params(level, rpm)
        p = {k: v for k, v in src.items() if k not in ("t_end", "n_mix_cycles", "n_rec_periods")}
        p.update(run_id=rid(level, rpm, "rep"), _binary=str(po.BIN["rep"]), n_mix_cycles=0,
                 t_checkpoint=src["n_mix_cycles"] * T, t_end=db["t_end"] - src["n_mix_cycles"] * T,
                 restart_continue=1, n_rec_periods=1, replay_loop=1, replay_scalar_level=level, replay_diag=0,
                 replay_project_flux=0)
        ph = rg.WORK / "prereg_placeholder.dump"
        ph.write_bytes(b"")
        queued = subprocess.run(["squeue", "-h", "-j", jid, "-o", "%i"], capture_output=True, text=True).stdout.strip()
        out[key] = po._submit(p, po.NT_A, "24:00:00", checkpoint=str(ph),
                              dependency=f"afterok:{jid}" if queued else None)
        staged = rg.scratch_dir(p["run_id"]) / "checkpoint.dump"
        staged.unlink()
        staged.symlink_to(rg.scratch_dir(src["run_id"]) / "replay_start.dump")
    print(out)


def report(level: int) -> None:
    for rpm in RPMS:
        run = rg.RUNS / rid(level, rpm, "rep")
        if not (run / "params.json").exists():
            continue
        r = json.loads((run / "results.json").read_text()) if (run / "results.json").exists() else pp.main(str(run))
        d = json.loads((ROOT / "runs" / DB[level].format(rpm) / "results.json").read_text())
        print(f"{rpm:5g} rpm  " + "  ".join(f"{q} DNS {d[q]:7.2f} replay {r[q]:7.2f} ({100 * (r[q] / d[q] - 1):+5.1f} %)"
                                          for q in ("kLa_10", "dtmix_0.50", "dtmix_0.75")))


if __name__ == "__main__":
    cmd, lev = sys.argv[1], int(sys.argv[2])
    if cmd == "diff":
        diff(lev, float(sys.argv[3]))
    elif cmd == "dns":
        dns(lev, [float(x) for x in sys.argv[3:]] or list(RPMS))
    else:
        {"replay": replay, "report": report}[cmd](lev)
