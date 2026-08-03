"""One-off submission script for the L9/L10 short-window fidelity test at 30.0 RPM.

Submits two cold-start, single-segment, non-chained SLURM jobs (fidelity 9 and
fidelity 10) truncated to t_end=8.0 (one period past the ramp). No checkpoint
restart is involved in either job.
"""
import json
import sys
from pathlib import Path
from uuid import uuid4

sys.path.insert(0, str(Path(__file__).parents[1]))

from scripts.simulate import submit_slurm

PROJECT_ROOT = Path(__file__).parents[1]
EXP_DIR = PROJECT_ROOT / "experiments" / "l9_l10_short_window_test_30rpm"
MPI_TEMPLATE = PROJECT_ROOT / "config" / "slurm_mpi_template.sh"

def main():
    jobs = {}
    for label, fname in [("f9_short", "params_f9_short.json"), ("f10_short", "params_f10_short.json")]:
        params = json.loads((EXP_DIR / fname).read_text())
        run_id = uuid4().hex[:8]
        params["run_id"] = run_id
        params["_experiment_dir"] = str(EXP_DIR)

        walltime = params.pop("_walltime")
        mem = params.pop("_mem")
        ntasks = params.pop("_ntasks")
        binary = params.pop("_binary")
        params["_binary"] = binary

        job_id = submit_slurm(
            params,
            project_root=PROJECT_ROOT,
            walltime=walltime,
            template=MPI_TEMPLATE,
            mem=mem,
            cpus=1,
            ntasks=ntasks,
        )
        jobs[label] = {"run_id": run_id, "job_id": job_id}
        print(f"{label}: run_id={run_id} job_id={job_id}")

    (EXP_DIR / "_submitted_jobs.json").write_text(json.dumps(jobs, indent=2))


# [PROJECT FIX, 2026-08-03] This filename ends in "_test.py", matching pytest's
# default `python_files = test_*.py *_test.py` collection glob -- with no
# __main__ guard, a plain `pytest` invocation from the repo root IMPORTED this
# module and submitted two real SLURM jobs as a side effect (discovered when
# it happened twice in one session; see diary.md 2026-08-03). The leading
# underscore was meant to signal "not a test," but pytest's glob doesn't care.
if __name__ == "__main__":
    main()
