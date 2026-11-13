"""One-off: postprocess the 9 stalled sweep_tau_theta7_l9 runs from raw .dat files.

No new simulation needed -- all 9 runs already reached t_end on scratch;
they just never got a results.json because the SLURM jobs were walltime-killed
before the job script's postprocess step ran.
"""
from pathlib import Path

from scripts import postprocess as pp

RUN_IDS = [
    "44133566", "0183ca21", "8994c04a", "b1b72f63", "32c996c6",
    "488db14b", "30fb2321", "77e91d1b", "44b0f74d",
]
SCRATCH_ROOT = Path("/oscar/scratch/eaguerov/mpi_runs")

for rid in RUN_IDS:
    run_dir = SCRATCH_ROOT / rid
    print(f"=== {rid} ===")
    try:
        res = pp.main(str(run_dir))
        print(f"  kLa_25={res.get('kLa_25')}  vor_mean={res.get('vor_mean')}")
    except Exception as e:
        print(f"  FAILED: {e!r}")
