"""L8 flow-recording runs: the kmix_l8 / T5 protocol and binary, plus uniform-grid field snapshots after release.

Purpose (pre-registered, diary 2026-10-10): same-run test of a transfer-operator mixing estimator (is the rpm
roughness of dtmix in the flow or in the tracer?), and an exact-rerun reproducibility check of dtmix.
Only frames_per_period, snap_start_cycle, snap_end_cycle differ from scripts/submit_l8_t5_t6.params(); the lean
binary has no VIDEOS output, so dt_video drives only the snapshot event and the dynamics are unchanged.
Snapshots: fields/snap_%06d.bin, planes u.x, u.y, omega, f, cs, c, c1, c2, c3, oxy on the NN x NN grid.

Usage: uv run python scripts/submit_l8_record.py --rpm R [R ...] [--theta 7] [--dry-run]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm  # noqa: E402
from scripts.submit_l8_t5_t6 import NTASKS, SPINUP, WALL, params  # noqa: E402

FRAMES, SNAP_FROM, SNAP_TO = 100, SPINUP - 2, SPINUP + 10  # 2 periods before release, 10 after


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rpm", nargs="+", type=float, required=True)
    ap.add_argument("--theta", type=float, default=7.0)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    for rpm in a.rpm:
        p = params(rpm, a.theta, f"rec_l8_rpm{rpm:g}_th{a.theta:g}")
        p.update(frames_per_period=FRAMES, snap_start_cycle=float(SNAP_FROM), snap_end_cycle=float(SNAP_TO))
        print(f"  {p['run_id']}: t_end {p['t_end']}, snapshots cycles {SNAP_FROM}-{SNAP_TO}, {FRAMES}/period")
        if a.dry_run:
            continue
        job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=WALL,
                           template=ROOT / "config" / "slurm_mpi_template.sh", cpus=1, ntasks=NTASKS, mem="4G")
        print(f"    job={job}")


if __name__ == "__main__":
    main()
