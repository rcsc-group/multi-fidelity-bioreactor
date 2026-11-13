"""Golden-run gate: prove that a change to src/ leaves every solver output byte-identical.

A refactor of the solver (moving events into headers, renaming, deleting dead branches) must not change
the physics. This tool builds the solver from two source trees, runs the same short cases with each, and
compares a SHA-256 of every output file.

    golden.py run   SRC_DIR WORK_DIR       build the variants from SRC_DIR, run all cases, write WORK_DIR/manifest.json
    golden.py compare A/manifest.json B/manifest.json      exit 1 if any output differs
    golden.py compile SRC_DIR WORK_DIR     check that the research switches (COMPILE_ONLY) still build

Run `golden.py run` on a compute node (tools/golden/golden.sbatch does this), never on a login node.

Cases (level 6, serial, OMP_NUM_THREADS=1 so that sums are taken in a fixed order):
    base      default build; release after 1 period, run 2 periods; ends with checkpoint.dump
    restart   default build; restart_continue from base's checkpoint.dump for half a period
    diag      -DDIAGNOSTICS=1  (pressure_diag.dat)
    lean      -DEXTRA_TRACERS=0
    record    -DRECORD_FLUX=1; records the face fluxes for 1 period after the release
    replay    -DREPLAY_FLUX=1; advects the scalars on record's fluxes for 1 period
Not covered: MPI builds (no MPI on this partition) and VIDEOS (needs an OpenGL context).

A file that differs between two runs of the SAME source (a wall-clock column, say) cannot be a gate. Run the
baseline twice and pass both manifests to `compare --drop-unstable BASE2`; such files are listed, not hidden.
"""
import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

QCC = "/oscar/data/dharri15/eaguerov/basilisk/src/qcc"
BASILISK = "/oscar/data/dharri15/eaguerov/basilisk/src"
CFLAGS = ["-O2", "-w", "-fopenmp", "-Wall", "-autolink"]  # as `make build`
LDFLAGS = [f"-L{BASILISK}/gl", "-lglutils", "-lfb_tiny", "-lGL", "-lm"]
BUILDS = {"default": [], "diag": ["-DDIAGNOSTICS=1"], "lean": ["-DEXTRA_TRACERS=0"],
          "record": ["-DRECORD_FLUX=1"], "replay": ["-DREPLAY_FLUX=1"]}
# Compile-only: switches the cases do not run. They must still build (a switch that no longer compiles is
# found the day someone needs it, which is too late). `golden.py compile SRC_DIR WORK_DIR`.
COMPILE_ONLY = {"videos": ["-DVIDEOS=1"], "crosslevel": ["-DCROSS_LEVEL_WARMSTART=1"],
                "crosslevel_old": ["-DCROSS_LEVEL_WARMSTART=1", "-DCROSS_LEVEL_PREREFINE_SOLID=0"],
                "kimflux": ["-DRECORD_FLUX=1", "-DREC_KIM_FLUX=1", "-DREC_KIM_MASK=1"],
                "prereg": ["-DRECORD_FLUX=1", "-DPREREG_DIAG=1"], "legacy_su": ["-DLEGACY_SU=1"],
                "checksum_rec": ["-DRECORD_FLUX=1", "-DREC_CHECKSUM=1"], "checksum_rep": ["-DREPLAY_FLUX=1", "-DREC_CHECKSUM=1"],
                "seed_f05": ["-DSEED_KIM=0"], "no_streaming": ["-DSTREAMING_VORTICITY=0"],
                "slow_ramp": ["-DN_RAMP_CYCLES=40"]}
SKIP = {"params.json", "stdout.txt", "stderr.txt"}  # inputs and logs; logs carry wall-clock times


def base_params(run_id: str) -> dict:
    z = [0.0, 0.0, 0.0]
    return {"run_id": run_id, "fidelity": 6, "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
            "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": z, "omega_b": 32.5 * 2 * math.pi / 60,
            "omega_h": 0.0, "amplitude_h": z, "phi_horizontal": z, "n_mix_cycles": 1}


def build(src: Path, work: Path, builds: dict = BUILDS) -> dict:
    env = {**os.environ, "BASILISK": BASILISK, "PATH": f"{BASILISK}:{os.environ['PATH']}"}
    bins = {}
    for name, flags in builds.items():
        d = work / "build" / name
        shutil.rmtree(d, ignore_errors=True)
        shutil.copytree(src, d)
        out = d / "BioReactor"
        r = subprocess.run([QCC, *CFLAGS, *flags, "BioReactor.c", "-o", str(out), *LDFLAGS],
                           cwd=d, env=env, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"build {name} failed:\n{r.stderr[-3000:]}")
        bins[name] = out
    return bins


def run(binary: Path, rd: Path, params: dict, restart: Path | None = None) -> None:
    shutil.rmtree(rd, ignore_errors=True)
    rd.mkdir(parents=True)
    (rd / "params.json").write_text(json.dumps(params, indent=1))
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    cmd = [str(binary), "params.json"] + ([str(restart)] if restart else [])
    with open(rd / "stdout.txt", "w") as o, open(rd / "stderr.txt", "w") as e:
        r = subprocess.run(cmd, cwd=rd, env=env, stdout=o, stderr=e)
    if r.returncode != 0:
        sys.exit(f"{rd.name}: exit {r.returncode}; see {rd}/stderr.txt")


def digest(rd: Path) -> dict:
    out = {}
    for f in sorted(rd.rglob("*")):
        if f.is_file() and f.name not in SKIP:
            out[str(f.relative_to(rd))] = hashlib.sha256(f.read_bytes()).hexdigest()
    return out


def cmd_run(src: Path, work: Path) -> None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from scripts.postprocess import _t_scales  # the solver's own period, as the campaign scripts use it
    work.mkdir(parents=True, exist_ok=True)
    bins = build(src, work)
    T = _t_scales(base_params("x"))[1]
    R = work / "runs"
    p = base_params("base") | {"t_end": 2 * T}
    run(bins["default"], R / "base", p)
    t_ck = 2 * T  # absolute end time of base (no checkpoint offset)
    run(bins["default"], R / "restart", base_params("restart") | {"t_checkpoint": t_ck, "restart_continue": 1,
                                                                  "n_mix_cycles": 0, "t_end": 0.5 * T},
        restart=R / "base" / "checkpoint.dump")
    run(bins["diag"], R / "diag", base_params("diag") | {"t_end": 1.5 * T})
    run(bins["lean"], R / "lean", base_params("lean") | {"t_end": 1.5 * T})
    run(bins["record"], R / "record", base_params("record") | {"t_end": 2 * T, "n_rec_periods": 1, "rec_dir": "."})
    run(bins["replay"], R / "replay", base_params("replay") | {"t_checkpoint": T, "restart_continue": 1,
                                                               "n_mix_cycles": 0, "t_end": T,
                                                               "rec_dir": str(R / "record")},
        restart=R / "record" / "replay_start.dump")
    man = {c.name: digest(c) for c in sorted(R.iterdir())}
    (work / "manifest.json").write_text(json.dumps(man, indent=1))
    print(f"{sum(len(v) for v in man.values())} files hashed in {len(man)} cases -> {work}/manifest.json")


def cmd_compare(a: Path, b: Path, unstable: Path | None) -> int:
    A, B = json.loads(a.read_text()), json.loads(b.read_text())
    drop = set()
    if unstable:
        U = json.loads(unstable.read_text())
        drop = {(c, f) for c in A for f in A[c] if U.get(c, {}).get(f) != A[c][f]}
        for c, f in sorted(drop):
            print(f"UNSTABLE (differs between two runs of the same source, not gated): {c}/{f}")
    bad = 0
    for c in sorted(set(A) | set(B)):
        for f in sorted(set(A.get(c, {})) | set(B.get(c, {}))):
            if (c, f) in drop:
                continue
            x, y = A.get(c, {}).get(f), B.get(c, {}).get(f)
            if x != y:
                bad += 1
                print(f"DIFF {c}/{f}: {'missing' if x is None else x[:12]} vs {'missing' if y is None else y[:12]}")
    n = sum(len(v) for v in A.values()) - len(drop)
    print(f"{'FAIL' if bad else 'PASS'}: {bad} of {n} gated files differ")
    return 1 if bad else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("src", type=Path); r.add_argument("work", type=Path)
    c = sub.add_parser("compare"); c.add_argument("a", type=Path); c.add_argument("b", type=Path)
    c.add_argument("--drop-unstable", type=Path, default=None)
    k = sub.add_parser("compile"); k.add_argument("src", type=Path); k.add_argument("work", type=Path)
    a = ap.parse_args()
    if a.cmd == "run":
        cmd_run(a.src.resolve(), a.work.resolve())
    elif a.cmd == "compile":
        build(a.src.resolve(), a.work.resolve(), COMPILE_ONLY)
        print(f"{len(COMPILE_ONLY)} compile-only variants build: {', '.join(COMPILE_ONLY)}")
    else:
        sys.exit(cmd_compare(a.a, a.b, a.drop_unstable))
