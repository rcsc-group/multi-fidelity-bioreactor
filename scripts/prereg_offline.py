"""Offline-vs-DNS test, our solver (pre-registered in diary.md, 2026-10-29). 20 rpm, theta 7, from rest, release at
cycle N_MIX, 3 recorded cycles.
  build             BioReactor-pr-rec (-DRECORD_FLUX), pr-kimrec (+ -DREC_KIM_FLUX=1), pr-rep (-DREPLAY_FLUX), MPI
  dns  LEVEL        A (record, NT_A ranks), Akim (Kim-style record, NT_A ranks), B (record, NT_B ranks: noise floor)
  replay LEVEL      afterok of the DNS jobs: T1 (A, no loop), H (A, hold 100), W (A, 1 cycle looped), S (Akim, no loop)
  report LEVEL      metrics M1-M3 of every run against A, written to docs/meetings/prereg/
Run ids: pr{LEVEL}_{tag}. Smoke test (gate G2): LEVEL 6 first.
Usage: uv run python scripts/prereg_offline.py build | dns 6 | replay 6 | report 6
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402
from scripts.simulate import submit_slurm  # noqa: E402

RPM, THETA = 20.0, 7.0
N_MIX, N_REC = 20, 3
NT_A, NT_B = 16, 8
EXCLUDE = "node2356"
REC_ROOT = Path("/oscar/data/dharri15/eaguerov/replay_records/prereg")
BIN = {"rec": rg.BIN_DIR / "BioReactor-pr4-rec", "kimrec": rg.BIN_DIR / "BioReactor-pr2-kimrec",
       "rep": rg.BIN_DIR / "BioReactor-pr4-rep",   # pr2/pr4: + per-step diag (-DPREREG_DIAG=1)
       "recseed": rg.BIN_DIR / "BioReactor-pr5-recseed", "kimmask": rg.BIN_DIR / "BioReactor-pr5-kimmask"}
FLAGS = {"rec": ["-DRECORD_FLUX", "-DEXTRA_TRACERS=0", "-DPREREG_DIAG=1"],
         "kimrec": ["-DRECORD_FLUX", "-DEXTRA_TRACERS=0", "-DREC_KIM_FLUX=1", "-DPREREG_DIAG=1"],
         "rep": ["-DREPLAY_FLUX", "-DEXTRA_TRACERS=0", "-DPREREG_DIAG=1"],
         "recseed": ["-DRECORD_FLUX", "-DEXTRA_TRACERS=0", "-DPREREG_DIAG=1", "-DSEED_KIM=1"],
         "kimmask": ["-DRECORD_FLUX", "-DEXTRA_TRACERS=0", "-DREC_KIM_FLUX=1", "-DREC_KIM_MASK=1", "-DPREREG_DIAG=1"]}
OUT = ROOT / "docs/meetings/prereg"
PFX = os.environ.get("PREREG_PFX", "pr")   # run-id prefix: "pr2" = rerun with the per-step diag


def rid(level: int, tag: str) -> str:
    return f"{PFX}{level}_{tag}"


def rec_dir(level: int, tag: str) -> Path:
    return REC_ROOT / rid(level, tag)


def base(level: int, run_id: str, binary: Path) -> dict:
    return rg.make_params(run_id, binary, fidelity=level, n_mix_cycles=N_MIX, post_cycles=N_REC + 0.05,
                          rpm=RPM, theta=THETA)


def _submit(params: dict, ntasks: int, walltime: str, checkpoint: str | None = None,
            dependency: str | None = None) -> str:
    for d in (rg.run_dir(params["run_id"]), rg.scratch_dir(params["run_id"])):
        if d.exists():
            shutil.rmtree(d)
    rg.RUNS.mkdir(parents=True, exist_ok=True)
    return submit_slurm(params, project_root=ROOT, runs_root=rg.RUNS, walltime=walltime,
                        template=ROOT / "config" / "slurm_mpi_template.sh", cpus=1, ntasks=ntasks, mem="4G",
                        checkpoint=checkpoint, dependency=dependency, exclude=EXCLUDE,
                        mpi_scratch_root=rg.MPI_SCRATCH)


def build(only: list[str] | None = None) -> None:
    for k, b in BIN.items():
        if only and k not in only:
            continue
        if b.exists():
            raise RuntimeError(f"{b} exists: never overwrite a binary; remove it by hand if it is yours")
        print(rg.build_binary(b.name.removeprefix("BioReactor-"), FLAGS[k], mpi=True))


def dns(level: int, only: list[str] | None = None) -> None:
    wall = "01:00:00" if level <= 6 else "12:00:00"
    f = rg.WORK / f"prereg{PFX[2:]}_dns_jobs_{level}.json"
    jobs = json.loads(f.read_text()) if only and f.exists() else {}
    for tag, b, nt in (("A", "rec", NT_A), ("Akim", "kimrec", NT_A), ("B", "rec", NT_B),
                       ("Aseed", "recseed", NT_A), ("Akimmask", "kimmask", NT_A)):
        if only and tag not in only:
            continue
        d = rec_dir(level, tag)
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
        p = base(level, rid(level, tag), BIN[b])
        p.update(n_rec_periods=N_REC, rec_dir=str(d))
        jobs[tag] = _submit(p, nt, wall)
    f.write_text(json.dumps(jobs))
    print(jobs)


def replay(level: int, only: list[str] | None = None) -> None:
    jobs = json.loads((rg.WORK / f"prereg{PFX[2:]}_dns_jobs_{level}.json").read_text())
    wall = "01:00:00" if level <= 6 else "06:00:00"
    out = {}
    for tag, src, extra in (("T1", "A", dict(n_rec_periods=N_REC, replay_loop=0)),
                            ("H", "A", dict(n_rec_periods=N_REC, replay_loop=0, replay_hold=100)),
                            ("W", "A", dict(n_rec_periods=1, replay_loop=1)),
                            ("S", "Akim", dict(n_rec_periods=N_REC, replay_loop=0)),
                            ("WS", "Akim", dict(n_rec_periods=1, replay_loop=1)),   # flux from u, cycle 1 looped
                            ("WSm", "Akimmask", dict(n_rec_periods=1, replay_loop=1))):   # + Kim's solid-face mask
        if only and tag not in only:
            continue
        sp = base(level, rid(level, f"rep{tag}"), BIN["rep"])
        _, T = pp._t_scales(sp)
        t_rec0 = N_MIX * T
        p = {k: v for k, v in sp.items() if k not in ("t_end", "n_mix_cycles")}
        # t_end of a replay counts from the checkpoint (as replay_lf.py: t_end = N_REPLAY_CYCLES * T)
        p.update(n_mix_cycles=0, t_checkpoint=t_rec0, t_end=(N_REC + 0.02) * T, restart_continue=1,
                 rec_dir=str(rec_dir(level, src)), replay_scalar_level=level, replay_diag=0, replay_project_flux=0,
                 **extra)
        # the record's replay_start.dump does not exist yet: stage a placeholder, then link it (as replay_lf.py)
        ph = rg.WORK / "prereg_placeholder.dump"
        ph.write_bytes(b"")
        queued = subprocess.run(["squeue", "-h", "-j", jobs[src], "-o", "%i"], capture_output=True, text=True).stdout.strip()
        if not queued and not (rg.scratch_dir(rid(level, src)) / "replay_start.dump").exists():
            raise RuntimeError(f"DNS {src} ({jobs[src]}) is not queued and left no replay_start.dump")
        jid = _submit(p, NT_A, wall, checkpoint=str(ph), dependency=f"afterok:{jobs[src]}" if queued else None)
        staged = rg.scratch_dir(p["run_id"]) / "checkpoint.dump"
        staged.unlink()
        staged.symlink_to(rg.scratch_dir(rid(level, src)) / "replay_start.dump")
        out[tag] = jid
    print(out)


def series(run_id: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Time (cycles from the release), tracer mass in the liquid, chi, liquid oxygen: as fig_kla_mix_compare.ours()."""
    run = rg.RUNS / run_id
    sp = json.loads((run / "params.json").read_text())
    _, T = pp._t_scales(sp)
    t, oxy, c, c2 = pp._joined_cols(run, "tr_oxy.dat", [2, 8, 9])
    tf, fl = pp._joined_cols(run, "vol_frac_interf.dat", [pp._COL_F_LIQ])
    fm = fl[(tf >= N_MIX * T - 1e-9 * T) & (tf <= (N_MIX + N_REC) * T)].mean()   # same window in every run
    m = t >= N_MIX * T - 1e-9 * T
    t, oxy, c, c2 = t[m], oxy[m], c[m], c2[m]
    k = int(np.argmax(c > 0))       # first row after the tracer is seeded
    t, oxy, c, c2 = t[k:], oxy[k:], c[k:], c2[k:]
    mean = c / fm
    var = c2 / fm - mean ** 2
    return (t - N_MIX * T) / T, c / c[0], 1 - var / var[0], oxy


def report(level: int) -> None:
    ta, ma, xa, oa = series(rid(level, "A"))
    rows, ts = [], {"t": ta, "A_m": ma, "A_chi": xa, "A_oxy": oa}
    for tag in ("B", "repT1", "repH", "repW", "repS", "repWS"):
        try:
            t, m, x, o = series(rid(level, tag))
        except (FileNotFoundError, OSError, ValueError) as e:
            rows.append(f"{tag}: missing ({e})")
            continue
        w = (ta > 0) & (ta <= min(t[-1], ta[-1], N_REC + 1e-9))
        mi, xi, oi = (np.interp(ta[w], t, v) for v in (m, x, o))
        m1 = np.max(np.abs(mi - ma[w]))
        m2 = np.max(np.abs(xi - xa[w]))
        w3 = ta[w] >= 0.25                       # M3 window [0.25, 3] cycles (diary amendment 2)
        m3 = np.max(np.abs(oi[w3] / oa[w][w3] - 1))
        rows.append(f"{tag:6s} M1 {m1:.3e}  M2 {m2:.3e}  M3 {m3:.3e}  (to {min(t[-1], ta[-1], N_REC):.3f} cycles)")
        for k, v in (("m", m), ("chi", x), ("oxy", o)):
            ts[f"{tag}_{k}"] = np.interp(ta, t, v, left=np.nan, right=np.nan)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"ours_L{level}_metrics.txt").write_text("\n".join(rows) + "\n")
    keys = list(ts)
    np.savetxt(OUT / f"ours_L{level}_series.csv", np.column_stack([ts[k] for k in keys]), delimiter=",",
               header=",".join(keys), comments="")
    print("\n".join(rows))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "build":
        build(sys.argv[2:] or None)                       # optional keys: rec kimrec rep
    elif cmd == "replay":
        replay(int(sys.argv[2]), sys.argv[3:] or None)     # optional tags: T1 H W S
    elif cmd == "dns":
        dns(int(sys.argv[2]), sys.argv[3:] or None)        # optional tags: A Akim B
    else:
        report(int(sys.argv[2]))
