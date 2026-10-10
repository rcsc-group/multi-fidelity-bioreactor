"""L8 flow-recording runs: the kmix_l8 / T5 protocol and binary, plus uniform-grid field snapshots after release.

Purpose (pre-registered, diary 2026-10-10): same-run test of a transfer-operator mixing estimator (is the rpm
roughness of dtmix in the flow or in the tracer?), and an exact-rerun reproducibility check of dtmix.
Only frames_per_period, snap_start_cycle, snap_end_cycle differ from scripts/submit_l8_t5_t6.params(); the lean
binary has no VIDEOS output, so dt_video drives only the snapshot event and the dynamics are unchanged.
Snapshots: fields/snap_%06d.bin, planes u.x, u.y, omega, f, cs, c, c1, c2, c3, oxy on the NN x NN grid.

Usage: uv run python scripts/submit_l8_record.py --rpm R [R ...] [--theta 7] [--level 8] [--frames 100]
       [--snap-from 78] [--snap-to 90] [--tag rec] [--end-cycle C] [--dry-run]
L9 uses 32 ranks and 30 h (fig9_l9_rpm32.5: 18.2 h at 32 ranks).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm  # noqa: E402
from scripts.submit_l8_t5_t6 import NTASKS, POST, SPINUP, WALL, params  # noqa: E402

FRAMES, SNAP_FROM, SNAP_TO = 100, SPINUP - 2, SPINUP + 10  # 2 periods before release, 10 after


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rpm", nargs="+", type=float, required=True)
    ap.add_argument("--theta", type=float, default=7.0)
    ap.add_argument("--level", type=int, default=8)
    ap.add_argument("--frames", type=int, default=FRAMES)
    ap.add_argument("--snap-from", type=float, default=SNAP_FROM)
    ap.add_argument("--snap-to", type=float, default=SNAP_TO)
    ap.add_argument("--tag", default="rec", help="run id prefix")
    ap.add_argument("--end-cycle", type=float, default=None, help="stop after this many cycles (default: full protocol)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    ntasks, wall = (NTASKS, WALL) if a.level == 8 else (32, "30:00:00")
    for rpm in a.rpm:
        p = params(rpm, a.theta, f"{a.tag}_l{a.level}_rpm{rpm:g}_th{a.theta:g}")
        p.update(fidelity=a.level, frames_per_period=a.frames, snap_start_cycle=float(a.snap_from),
                 snap_end_cycle=float(a.snap_to))
        if a.end_cycle is not None:
            p["t_end"] = round(p["t_end"] / (SPINUP + POST) * a.end_cycle, 4)
        print(f"  {p['run_id']}: t_end {p['t_end']}, snapshots cycles {a.snap_from:g}-{a.snap_to:g}, "
              f"{a.frames}/period, {ntasks} ranks, {wall}")
        if a.dry_run:
            continue
        job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=wall,
                           template=ROOT / "config" / "slurm_mpi_template.sh", cpus=1, ntasks=ntasks, mem="4G")
        print(f"    job={job}")


if __name__ == "__main__":
    main()
