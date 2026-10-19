"""Helpers for the record/replay gate (diary 2026-10-16: F0 no regression, F1 machinery exact, F2 looping).

Build, submit and compare the three binaries of the gate:
  default       : unchanged coupled solver (F0 compares it before/after the change).
  -DRECORD_FLUX : coupled solver that also writes the face fluxes uf and a start checkpoint
                  (replay_start.dump, at the end of the step that lands on t_rec0; the first record is the next step).
  -DREPLAY_FLUX : no Navier-Stokes; advects f and the scalars on the recorded fluxes.

Record format (one file per rank, <rec_dir>/flux_rec.<pid>.bin; written by src/BioReactor.c):
  file header : char magic[8] = "FLXREC01", int32 pid, int32 npe
  per step    : int32 i, float64 t, float64 dt, int64 n, then n pairs (uf.x, uf.y) float64
                (the left and bottom face of each leaf cell, in the traversal order of foreach_cell)
All integers/doubles are little endian, no padding.
"""
from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402
from scripts.simulate import submit_slurm  # noqa: E402

BIN_DIR = Path("/oscar/data/dharri15/eaguerov/bin/replay-gate")
WORK = Path("/oscar/scratch/eaguerov/replay_gate")
RUNS = WORK / "runs"            # canonical run dirs (results.json lands here)
MPI_SCRATCH = WORK / "mpi_runs"  # where the job runs and writes the record
BASILISK = Path("/oscar/data/dharri15/eaguerov/basilisk/src")
MAGIC = b"FLXREC01"
FILE_HEADER_BYTES = 16
STEP_HEADER = np.dtype([("i", "<i4"), ("t", "<f8"), ("dt", "<f8"), ("n", "<i8")])  # 28 bytes, packed


# ── build ────────────────────────────────────────────────────────────────────

def source_mtime(src: Path = ROOT / "src") -> float:
    return max(p.stat().st_mtime for p in src.iterdir() if p.suffix in (".c", ".h"))


def binary_is_current(binary: Path, src: Path = ROOT / "src") -> bool:
    """True if the binary is newer than every source file (stale binaries cost two days once)."""
    return binary.exists() and binary.stat().st_mtime > source_mtime(src)


SBUILD = Path("/oscar/scratch/eaguerov/replay_gate/sbuild2.sh")


def build_binary(tag: str, flags: list[str] | None = None, src: Path = ROOT / "src", mpi: bool = False) -> Path:
    """Compile in a private copy of src on a COMPUTE node through sbatch (default account; never the
    interactive allocation). Same command line as `make build` (or the MPI one with mpi=True)."""
    out = BIN_DIR / f"BioReactor-{tag}"
    work = WORK / f"build_{tag}"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    for p in src.iterdir():
        if p.is_file():
            shutil.copy2(p, work / p.name)
    BIN_DIR.mkdir(parents=True, exist_ok=True)
    if out.exists():
        out.unlink()
    env = {**os.environ, "BUILD_MPI": "1" if mpi else "0"}
    r = subprocess.run(["sbatch", "--wait", "--export=ALL", str(SBUILD), str(work), str(out), *(flags or [])],
                       env=env, capture_output=True, text=True)
    if r.returncode != 0 or not out.exists():
        raise RuntimeError(f"build of {tag} failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    return out


# ── parameters ───────────────────────────────────────────────────────────────

def make_params(run_id: str, binary: Path, *, fidelity: int, n_mix_cycles: int, post_cycles: float,
                rpm: float = 32.5, theta: float = 7.0, **extra) -> dict:
    """The params() of scripts/submit_l8_t5_t6.py with a configurable level and length."""
    zeros = [0.0, 0.0, 0.0]
    omega_b = rpm * 2 * math.pi / 60.0
    p = {"run_id": run_id, "fidelity": fidelity, "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0},
         "fill_level": 0.5, "n_harmonics": 1, "theta_max": [theta, 0.0, 0.0], "phi_angular": zeros,
         "omega_b": omega_b, "omega_h": 0.0, "amplitude_h": zeros, "phi_horizontal": zeros,
         "n_mix_cycles": n_mix_cycles, "_binary": str(binary)}
    _, t_per = pp._t_scales(p)
    p["t_end"] = (n_mix_cycles + post_cycles) * t_per
    p.update(extra)
    return p


def t_period(params: dict) -> float:
    return pp._t_scales(params)[1]


# ── jobs ─────────────────────────────────────────────────────────────────────

def run_dir(run_id: str) -> Path:
    return RUNS / run_id


def scratch_dir(run_id: str) -> Path:
    return MPI_SCRATCH / run_id


def submit(params: dict, checkpoint: Path | None = None, walltime: str = "02:00:00", mem: str = "6G",
           src: Path = ROOT / "src", ntasks: int = 1, exclude: str | None = None) -> str:
    """One rank, default (priority) account: submit_slurm is called without account/qos on purpose.
    Refuses to submit a binary that is older than the sources it was built from."""
    binary = Path(params["_binary"])
    if not binary_is_current(binary, src):
        raise RuntimeError(f"{binary} is older than {src}: rebuild before submitting")
    for d in (run_dir(params["run_id"]), scratch_dir(params["run_id"])):
        if d.exists():
            shutil.rmtree(d)
    RUNS.mkdir(parents=True, exist_ok=True)
    return submit_slurm(params, project_root=ROOT, runs_root=RUNS, walltime=walltime,
                        template=ROOT / "config" / "slurm_mpi_template.sh", cpus=1, ntasks=ntasks, mem=mem,
                        checkpoint=str(checkpoint) if checkpoint else None, mpi_scratch_root=MPI_SCRATCH, exclude=exclude)


def wait(job: str, poll: float = 20.0, timeout: float = 6 * 3600) -> str:
    """Block until the job leaves the queue; return its final sacct state."""
    t0 = time.monotonic()
    while time.monotonic() - t0 < timeout:
        out = subprocess.run(["squeue", "-h", "-j", job, "-o", "%T"], capture_output=True, text=True).stdout.strip()
        if not out:
            break
        time.sleep(poll)
    else:
        raise TimeoutError(job)
    for _ in range(30):
        st = subprocess.run(["sacct", "-j", job, "-X", "-n", "-o", "State,Account"], capture_output=True,
                            text=True).stdout.split()
        if st and st[0] not in ("RUNNING", "PENDING"):
            if "mbessa" in " ".join(st[1:]):
                raise RuntimeError(f"job {job} ran under a forbidden account: {st}")
            return st[0]
        time.sleep(5)
    raise TimeoutError(f"sacct never reported a final state for {job}")


# ── comparison ───────────────────────────────────────────────────────────────

def load_tr_oxy(path: Path) -> np.ndarray:
    return pp._load_dat(path)


def rel_diff_by_time(a: np.ndarray, b: np.ndarray, cols: list[int], t_tol: float = 1e-6) -> dict:
    """Relative difference of columns `cols` at every time present in both tables (matched on column 1 = t).

    rel = |a-b| / max(|a|, |b|, tiny); rows whose both values are exactly 0 are equal (rel 0).
    Returns {"n_common", "n_a", "n_b", "max": {col: value}, "t_common": array}.
    """
    ta, tb = a[:, 1], b[:, 1]
    idx_b = np.searchsorted(tb, ta)
    pairs = []
    for ia, t in enumerate(ta):
        for ib in (idx_b[ia] - 1, idx_b[ia]):
            if 0 <= ib < len(tb) and abs(tb[ib] - t) <= t_tol * max(1.0, abs(t)):
                pairs.append((ia, ib))
                break
    ia = np.array([p[0] for p in pairs], dtype=int)
    ib = np.array([p[1] for p in pairs], dtype=int)
    out = {}
    for c in cols:
        x, y = a[ia, c], b[ib, c]
        den = np.maximum(np.maximum(np.abs(x), np.abs(y)), 1e-300)
        out[c] = float(np.max(np.abs(x - y) / den)) if len(ia) else math.nan
    return {"n_common": len(ia), "n_a": len(ta), "n_b": len(tb), "max": out, "t_common": ta[ia]}


# ── record reader ────────────────────────────────────────────────────────────

def read_record(rec_dir: Path, pid: int = 0) -> tuple[np.ndarray, np.memmap | None, int]:
    """Memory-map a flux record. Returns (headers, data, n) with data[k] = (4n,) floats of step k, or
    (headers, None, 0) if the record is empty."""
    path = Path(rec_dir) / f"flux_rec.{pid}.bin"
    with open(path, "rb") as fh:
        if fh.read(8) != MAGIC:
            raise ValueError(f"{path}: bad magic")
    raw = np.memmap(path, dtype=np.uint8, mode="r", offset=FILE_HEADER_BYTES)
    if raw.size == 0:
        return np.empty(0, STEP_HEADER), None, 0
    n = int(np.frombuffer(raw[20:28].tobytes(), "<i8")[0])
    rec = STEP_HEADER.itemsize + 32 * n
    nrec = raw.size // rec
    if nrec * rec != raw.size:
        raise ValueError(f"{path}: size {raw.size} is not a multiple of the record size {rec}")
    mm = np.memmap(path, dtype=np.uint8, mode="r", offset=FILE_HEADER_BYTES, shape=(nrec, rec))
    heads = np.empty(nrec, STEP_HEADER)
    for k in range(nrec):
        heads[k] = np.frombuffer(mm[k, :STEP_HEADER.itemsize].tobytes(), STEP_HEADER)[0]
    return heads, mm, n


def step_flux(mm: np.memmap, k: int, n: int) -> np.ndarray:
    return np.frombuffer(mm[k, STEP_HEADER.itemsize:].tobytes(), "<f8")


def flow_defect(rec_dir: Path, t_per: float, n_periods: int) -> dict:
    """max |uf(period k) - uf(period 1)| / max|uf| at equal phase, for k = 2..n_periods.

    The timestep sequence differs from period to period, so uf of period k is interpolated linearly in
    time to the phases (t_j - t_0) of the steps j of period 1. Returns the worst value over all k
    and faces, plus the value per k.
    """
    heads, mm, n = read_record(rec_dir)
    t = heads["t"]
    t0 = t[0]
    j1 = np.nonzero(t < t0 + t_per * (1 - 1e-9))[0]
    umax = max(np.max(np.abs(step_flux(mm, k, n))) for k in range(len(heads)))
    per_k = []
    for k in range(2, n_periods + 1):
        worst = 0.0
        for j in j1:
            tt = t[j] + (k - 1) * t_per
            hi = int(np.searchsorted(t, tt))
            if hi <= 0 or hi >= len(t):
                break
            lo = hi - 1
            w = (tt - t[lo]) / (t[hi] - t[lo])
            u = (1 - w) * step_flux(mm, lo, n) + w * step_flux(mm, hi, n)
            worst = max(worst, float(np.max(np.abs(u - step_flux(mm, j, n)))))
        per_k.append(worst / umax)
    return {"per_period": per_k, "max": max(per_k) if per_k else math.nan, "umax": float(umax)}


def cstar_series(run: Path) -> tuple[np.ndarray, np.ndarray]:
    return pp._compute_c_star(run)


def load_results(run: Path) -> dict:
    return json.loads((run / "results.json").read_text())
