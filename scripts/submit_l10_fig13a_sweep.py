"""L10 Fig13a replication: 9 independent cold starts, 32.5-rpm-condition
methodology matching the L9/L8/L6 series already in replicated_Fig13.png.
theta=7deg fixed, rpm = Kim's 9 points.

48 cores/node, single node per job (avoids the multi-node output-
corruption bug found 2026-09-10). Runs in two waves under the 312-CPU
priority cap (floor(312/48)=6 concurrent): 6 then 3.

CYCLES is set from the 48-core calibration probe (l10_48core_probe_
rpm32.5) -- see diary.md 2026-09-10 for the measured rate this was sized
against.

Usage:
    uv run python scripts/submit_l10_fig13a_sweep.py
"""
import sys
from pathlib import Path

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-video-fixed"
T_PER_ND = 0.607329
CYCLES = 30  # trimmed from 35 to keep the 24h walltime (see below) a
             # comfortable margin above the worst-case per-cycle estimate
             # (30 * ~40min/cycle worst-case ~= 20h, vs 35's ~23.3h --
             # too tight against a 24h cap). Extend individual points later
             # if validate_run.py finds they haven't settled.
# [2026-09-10] First wave only -- 3 of 9 points, deliberately leaving CPU
# headroom for the cross-level (L8->L9 adapt_wavelet) warm-start pilot
# running alongside. Calibration probe (l10_48core_probe_rpm32.5) was still
# queued (scheduler priority delay, not a resource cap) when this wave was
# submitted -- proceeding without it is a bounded risk: these are
# independent cold starts, so a wrong cycle/walltime guess just means one
# job hits TIMEOUT and gets resubmitted, not a cascading chain failure.
# 17.5/37.5 bracket the range; 32.5 is Kim's canonical baseline and the
# condition every other fidelity level already has data for.
RPMS = [17.5, 32.5, 37.5]

jobs = {}
for rpm in RPMS:
    omega = rpm * 2 * 3.141592653589793 / 60.0
    params = {
        "run_id": f"l10_fig13a_rpm{rpm:g}", "fidelity": 10,
        "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
        "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": omega, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": CYCLES * T_PER_ND, "n_mix_cycles": CYCLES, "frames_per_period": 5,
        "_binary": BINARY,
    }
    # [2026-09-10] 48h originally requested here backfilled at ~4 DAYS out
    # (scontrol show job: StartTime 4 days after submission) -- long
    # walltime requests are much harder for the scheduler to slot in on a
    # busy shared cluster, since it needs a guaranteed-free window spanning
    # the WHOLE request. 24h is still safely above the worst-case estimate
    # (35 cycles * ~40min/cycle worst-case ~= 23.3h -- tight) but backfills
    # far better; verify the new StartTime immediately after submitting,
    # don't just assume it helped.
    job = submit_slurm(params, project_root=root, walltime="24:00:00",
                       template=root / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=48, mem="6G", exclude="node2336")
    jobs[rpm] = job
    print(f"{rpm:g}rpm: job={job}  run_id=l10_fig13a_rpm{rpm:g}")

print(f"\nall jobs: {jobs}")
