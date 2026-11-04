"""Submit the edr_95 backfill restarts (lead, 2026-11-04): 3 cycles from the last checkpoint of each run in
state/edr_backfill_list.csv, binary BioReactor-mpi-lean-2ef5b07, new run id <source>_e95 (source folders untouched).

Each sbatch command is printed in full before it is sent. Checks per job: partition batch, account default,
--cpus-per-task=1, --ntasks=<ranks of the level>, --mem-per-cpu=4G (MPI template), --exclude=node2356, no mbessa.
Usage: uv run python scripts/campaign2d/submit_edr95_backfill.py [--go]     (default: print only)
"""
import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts import simulate  # noqa: E402
from scripts.submit_continue import build_params  # noqa: E402
from scripts.postprocess import _t_scales  # noqa: E402

BINARY = "/oscar/data/dharri15/eaguerov/bin/BioReactor-mpi-lean-2ef5b07"
LIST = ROOT / "scripts" / "campaign2d" / "state" / "edr_backfill_list.csv"
OUT = ROOT / "scripts" / "campaign2d" / "state" / "jobs_edr95.jsonl"
MPI_SCRATCH = Path("/oscar/scratch/eaguerov/mpi_runs")
RANKS = {8: 16, 9: 32, 10: 32}                       # plan 2 (STAGES)
WALL = {8: "04:00:00", 9: "08:00:00", 10: "12:00:00"}  # 3 cycles; L10 at 42-168 min per cycle, 32 ranks
N_CYC = 3
FORBIDDEN = "mbessa-condo"


def printing_run(cmd, **kw):
    """Print the full sbatch command, refuse forbidden or wrong flags, then run it."""
    line = " ".join(cmd)
    print("SBATCH:", line)
    assert FORBIDDEN not in line, "mbessa-condo is forbidden"
    assert "--cpus-per-task=1" in line and "--exclude=node2356" in line and "--partition=batch" in line
    assert "--account=default" in line and "--mem-per-cpu=4G" in line, "flag check failed"
    return _ORIG_RUN(cmd, **kw)


_ORIG_RUN = subprocess.run


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--go", action="store_true")
    a = ap.parse_args()
    rows = list(csv.DictReader(LIST.open()))
    print(f"{len(rows)} restarts, binary {BINARY}")
    if a.go:
        simulate.subprocess.run = printing_run
    for r in rows:
        lev = int(r["level"])
        src = r["run"]
        sp = ROOT / "runs" / src
        p = build_params(json.loads((sp / "params.json").read_text()), src, "e95", _t_last(sp), N_CYC, 0.0, BINARY)
        rid = p["run_id"]
        if not rid.endswith("_e95") or rid == src:
            raise SystemExit(f"bad run id {rid}")
        if (ROOT / "runs" / rid).exists():
            print(f"skip {rid}: exists")
            continue
        stage = MPI_SCRATCH / rid
        if not a.go:
            print(f"  (dry run: {rid} would stage to {stage} and sbatch; nothing written)")
            continue
        if a.go:
            stage.mkdir(parents=True, exist_ok=True)
            if (sp / "checkpoint.phase").exists():
                shutil.copy2(sp / "checkpoint.phase", stage / "checkpoint.phase")
        print(f"{rid}: level {lev}, ranks {RANKS[lev]}, wall {WALL[lev]}, t_checkpoint {p['t_checkpoint']:.4f}")
        job = simulate.submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=WALL[lev],
                                    template=ROOT / "config" / "slurm_mpi_template.sh",
                                    checkpoint=str((sp / "checkpoint.dump").resolve()), cpus=1, ntasks=RANKS[lev],
                                    mem="4G", exclude="node2356", partition="batch", account="default")
        if a.go:
            with OUT.open("a") as f:
                f.write(json.dumps(dict(job=job, run=rid, source=src, level=lev, ranks=RANKS[lev])) + "\n")
            print(f"  job {job}")


def _t_last(sp):
    import numpy as np
    return float(np.loadtxt(sp / "shear_stress.dat", skiprows=1)[-1, 1])


if __name__ == "__main__":
    main()
