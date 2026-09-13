"""L10 Fig13a resubmit at 32 ranks instead of 48. The original 48-rank
attempts (l10_fig13a_rpm17.5/32.5/37.5, 2026-09-11) all TIMED OUT at 24h
with no checkpoint saved (BioReactor only dumps once, at the planned end)
-- that compute is unrecoverable, these are genuine fresh cold starts.

2026-09-13: found that 48-core (and even 64-core) single-node requests on
`batch` compete for a nearly-nonexistent pool of large-core nodes (checked
directly via sinfo: at the time, max free capacity on any single node was
~33-38 cores), while 32-core requests fit into much more commonly
available fragments and schedule far faster in practice. Real 32-rank L10
throughput measured today for 32.5rpm (Fig 8 late-time-check probe, job
6314896): 98.04 min/cycle.

No real 32-rank measurement exists yet for 17.5/32.5's siblings 17.5rpm or
37.5rpm -- rather than extrapolate and risk another blind-guess timeout,
this submits 32.5rpm's FULL 30-cycle target now (real data, safely
margined) and short 5-cycle CALIBRATION runs for 17.5/37.5rpm to get real
rates before committing a large walltime to them.

Usage:
    uv run python scripts/submit_l10_fig13a_32rank.py
"""
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-video-fixed"
T_PER_ND = 0.607329

# (rpm, n_mix_cycles, walltime, note)
PLAN = [
    (32.5, 30, "64:00:00", "full target -- real measured rate x1.3 margin"),
    (17.5, 5,  "24:00:00", "CALIBRATION only -- no real 32-rank rate yet, "
                            "sized from a heavily-margined scaled estimate"),
    (37.5, 5,  "14:00:00", "CALIBRATION only -- no real 32-rank rate yet, "
                            "sized from a heavily-margined scaled estimate"),
]

for rpm, cycles, walltime, note in PLAN:
    omega = rpm * 2 * 3.141592653589793 / 60.0
    run_id = f"l10_fig13a_rpm{rpm:g}_32rank" if cycles == 30 else f"l10_fig13a_rpm{rpm:g}_32rank_calib"
    params = {
        "run_id": run_id, "fidelity": 10,
        "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
        "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": omega, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": cycles * T_PER_ND, "n_mix_cycles": cycles, "frames_per_period": 5,
        "_binary": BINARY,
    }
    job = submit_slurm(params, project_root=root, runs_root=root / "runs",
                        walltime=walltime,
                        template=root / "config" / "slurm_mpi_template.sh",
                        cpus=1, ntasks=32, mem="4G")
    print(f"rpm={rpm:g}  run_id={run_id}  cycles={cycles}  walltime={walltime}  "
          f"job={job}  ({note})")
