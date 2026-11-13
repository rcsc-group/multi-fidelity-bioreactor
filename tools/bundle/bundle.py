"""Feature bundle: build every solver variant once, run the L6 cases of cases.py, write bundle.json.

    bundle.py run SRC_DIR WORK_DIR [--only CASE ...]     (on a compute node: tools/bundle/bundle.sbatch)

SRC_DIR is a solver source folder (src/). The tests in tests/bundle/ read WORK_DIR/bundle.json; they never
run anything. A case that fails is recorded with its error, its children are recorded as "not run", and the
other cases still run. The MPI case is run by its own job (MPI is not available on every partition).
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.bundle import cases as bc  # noqa: E402
from tools.golden import golden  # noqa: E402

MAXRUNTIME_RESERVE_S = 300   # src/BioReactor.c: the early dump fires this many seconds before the limit


def src_hash(src: Path) -> str:
    """SHA-256 over the solver sources (name + content, sorted): identity of a source tree without git."""
    h = hashlib.sha256()
    for f in sorted(list(src.glob("*.c")) + list(src.glob("*.h"))):
        h.update(f.name.encode() + b"\0" + f.read_bytes() + b"\0")
    return h.hexdigest()


def git_commit(src: Path):
    r = subprocess.run(["git", "-C", str(src), "rev-parse", "HEAD"], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else None


def dump_time(dump: Path) -> float:
    """Absolute time stored in the first 8 bytes of a Basilisk dump (what the solver reads on restart)."""
    import struct
    return struct.unpack("<d", dump.read_bytes()[:8])[0]


def resume_params(p: dict, dump: Path, t_end_abs: float) -> dict:
    """A continuation from `dump` that runs to the absolute time t_end_abs (t_end is a duration on restart)."""
    t_ck = dump_time(dump)
    return p | {"t_checkpoint": t_ck, "t_end": t_end_abs - t_ck}


def default_run_case(case, binary: Path, rd: Path, params: dict, restart, env: dict) -> None:
    rd.mkdir(parents=True, exist_ok=True)
    (rd / "params.json").write_text(json.dumps(params, indent=1))
    cmd = [str(binary), "params.json"] + ([str(restart)] if restart else [])
    with open(rd / "stdout.txt", "w") as o, open(rd / "stderr.txt", "w") as e:
        r = subprocess.run(cmd, cwd=rd, env={**os.environ, "OMP_NUM_THREADS": "1", **env}, stdout=o, stderr=e)
    if r.returncode != 0:
        tail = (rd / "stderr.txt").read_text()[-500:]
        raise RuntimeError(f"exit {r.returncode}: {tail}")


def default_build(src: Path, work: Path, builds: dict) -> dict:
    """One golden.build call per variant, so one failing variant does not stop the others."""
    out = {}
    for name, flags in builds.items():
        try:
            out.update(golden.build(src, work, {name: flags}))
        except SystemExit as e:
            out[name] = RuntimeError(f"build {name} failed: {str(e)[-500:]}")
    return out


def default_qoi(rd: Path) -> dict:
    from scripts import postprocess
    res = postprocess.main(str(rd))
    return {k: float(v) for k, v in res.items() if isinstance(v, (int, float)) and not isinstance(v, bool)}


def _with_parents(only, cases) -> dict:
    keep = set()
    for n in only:
        while n is not None:
            keep.add(n)
            n = cases[n].parent
    if "kill" in keep:
        keep.add("cold")   # kill's wall-clock limit comes from cold's wall time
    return {n: c for n, c in cases.items() if n in keep}


def run_bundle(src: Path, work: Path, only=None, cases=None, run_case=default_run_case,
               build_fn=default_build, qoi_fn=default_qoi, fidelity=None) -> dict:
    cases = cases or bc.CASES
    if only:
        cases = _with_parents(only, cases)
    work.mkdir(parents=True, exist_ok=True)
    runs = work / "runs"
    serial = {c.build for c in cases.values() if c.build != "mpi"}
    bins = build_fn(src, work, {b: bc.BUILDS[b] for b in sorted(serial)})
    out = {"src_hash": src_hash(src), "git_commit": git_commit(src), "builds": {}, "cases": {}}
    for b, path in bins.items():
        ok = not isinstance(path, Exception)
        sha = hashlib.sha256(Path(path).read_bytes()).hexdigest() if ok and Path(path).exists() else None
        out["builds"][b] = {"flags": bc.BUILDS[b], "sha256": sha, "error": "" if ok else str(path)}

    for name in bc.order(cases):
        c = cases[name]
        rd = runs / name
        rec = {"status": "not run", "error": "", "wall_s": 0.0, "dir": str(rd), "files": {}, "qoi": {}}
        out["cases"][name] = rec
        if c.build == "mpi":
            rec["error"] = "MPI case: runs in its own MPI job (tools/bundle/bundle_mpi.sbatch)"
            continue
        if c.parent and out["cases"][c.parent]["status"] != "ok":
            rec["error"] = f"parent {c.parent} failed"
            continue
        if isinstance(bins.get(c.build), Exception):
            rec["error"] = str(bins[c.build])
            continue
        p = bc.params(name, runs)
        if fidelity is not None:   # level check only (plan Task 2 step 6); the bundle itself is L6
            p["fidelity"] = fidelity
        env = dict(c.env)
        if env.get("BIOREACTOR_MAXRUNTIME_S") == "0":   # dump at about half of cold's wall time
            cw = out["cases"].get("cold", {})
            if cw.get("status") != "ok":
                rec["error"] = "needs cold's wall time (cold not ok)"
                continue
            env["BIOREACTOR_MAXRUNTIME_S"] = str(int(MAXRUNTIME_RESERVE_S + 0.5 * cw["wall_s"]))
        restart = runs / c.parent / c.restart if c.restart else None
        if name == "kill_resume":
            p = resume_params(p, restart, t_end_abs=3 * bc.period(bc.RPM_A))
        t0 = time.time()
        try:
            run_case(c, bins[c.build], rd, p, restart, env)
            rec["status"] = "ok"
        except Exception as e:  # noqa: BLE001 - recorded, not hidden
            rec["status"], rec["error"] = "failed", str(e)
        rec["wall_s"] = time.time() - t0
        if rd.exists():
            rec["files"] = golden.digest(rd)
        if rec["status"] == "ok":
            try:
                rec["qoi"] = qoi_fn(rd)
            except Exception as e:  # noqa: BLE001
                rec["qoi_error"] = f"{type(e).__name__}: {e}"
    (work / "bundle.json").write_text(json.dumps(out, indent=1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("src", type=Path)
    ap.add_argument("work", type=Path)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--fidelity", type=int, help="override the level of every case (level check only)")
    a = ap.parse_args()
    out = run_bundle(a.src.resolve(), a.work.resolve(), a.only, fidelity=a.fidelity)
    for n, c in out["cases"].items():
        print(f"{n:12s} {c['status']:8s} {c['wall_s']:7.1f} s  {c['error'][:100]}")
    sys.exit(0 if all(c["status"] == "ok" or "MPI" in c["error"] for c in out["cases"].values()) else 1)


if __name__ == "__main__":
    main()
