"""One-off: postprocess the 3 stalled L10 seg0 runs from already-complete raw data.

These jobs wrote their final checkpoint.dump hours ago and have produced zero
new output since (confirmed via file mtimes), while still burning CPU -- the
binary appears to hang past t_dump_checkpoint instead of exiting.  The raw
.dat trajectories are already complete and stable, so postprocessing now
(rather than waiting for a SLURM walltime kill) is safe.
"""
from pathlib import Path

from scripts import postprocess as pp

RUN_IDS = ["61d53bc5", "f02f61c2", "25cc70a4"]
SCRATCH_ROOT = Path("/oscar/scratch/eaguerov/mpi_runs")

for rid in RUN_IDS:
    run_dir = SCRATCH_ROOT / rid
    print(f"=== {rid} ===")
    try:
        res = pp.main(str(run_dir))
        print(f"  tau_98_qss={res.get('tau_98_qss')}  tau_100_max={res.get('tau_100_max')}  kLa_25={res.get('kLa_25')}")
    except Exception as e:
        print(f"  FAILED: {e!r}")
