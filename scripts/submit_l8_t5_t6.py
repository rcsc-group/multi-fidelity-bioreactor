"""L8 runs for tests T5 (dtmix roughness in rpm) and T6 (rpm x angle design), pre-registered in diary 2026-10-08.

Same protocol and binary as the kmix_l8 ladder (scripts/submit_kla_mix_ladder.py): cold start, 80-cycle
spin-up, tracer + oxygen released at cycle 80, then 150 cycles. Binary f1c11e0 (lean), copied to persistent
storage. Only rpm and theta change.
  T5: theta 7, rpm 20.625..31.875 step 0.625 (not already on the 2.5 grid) + rpm 22.4/22.6/29.9/30.1.
  T6: theta in {2, 4, 5.5, 9} x rpm 15..37.5 step 2.5.
Run ids: t5_l8_rpm<r>_th7, t6_l8_rpm<r>_th<theta>.

Usage: uv run python scripts/submit_l8_t5_t6.py --test t5|t6 [--dry-run] [--check RUN]
  --check RUN: print the params this script would write for RUN's (rpm, theta) and diff them with RUN's params.json.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm  # noqa: E402

BINARY = "/oscar/data/dharri15/eaguerov/bin/BioReactor-mpi-lean-f1c11e0"
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
SPINUP, POST, LEVEL, NTASKS, WALL = 80, 150, 8, 16, "08:00:00"
GRID = [15 + 2.5 * i for i in range(10)]


def params(rpm: float, theta: float, run_id: str) -> dict:
    zeros = [0.0, 0.0, 0.0]
    omega_b = rpm * 2 * math.pi / 60.0
    T_per = 2 * math.pi / omega_b
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    T_nd = T_per / (L / (V / (H * 0.5) / T_per))
    return {"run_id": run_id, "fidelity": LEVEL, "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
            "theta_max": [theta, 0.0, 0.0], "phi_angular": zeros, "omega_b": omega_b, "omega_h": 0.0,
            "amplitude_h": zeros, "phi_horizontal": zeros, "n_mix_cycles": SPINUP,
            "t_end": round((SPINUP + POST) * T_nd, 4), "_binary": BINARY}


def designs(test: str) -> list[tuple[float, float]]:
    if test == "t5":
        dense = [20 + 0.625 * i for i in range(21)]
        return [(r, 7.0) for r in dense if min(abs(r - g) for g in GRID) > 1e-9] + \
               [(r, 7.0) for r in (22.4, 22.6, 29.9, 30.1)]
    return [(r, th) for th in (2.0, 4.0, 5.5, 9.0) for r in GRID]


def check(run: str) -> None:
    old = json.load(open(ROOT / "runs" / run / "params.json"))
    rpm = old["omega_b"] * 60 / (2 * math.pi)
    new = params(rpm, old["theta_max"][0], old["run_id"])
    for k in sorted(set(old) | set(new)):
        if old.get(k) != new.get(k):
            print(f"  DIFF {k}: existing={old.get(k)!r} script={new.get(k)!r}")
    print(f"checked {run} (rpm {rpm:g}, theta {old['theta_max'][0]:g})")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", choices=("t5", "t6"))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--check")
    a = ap.parse_args()
    if a.check:
        check(a.check)
        return
    for rpm, th in designs(a.test):
        rid = f"{a.test}_l8_rpm{rpm:g}_th{th:g}"
        p = params(rpm, th, rid)
        print(f"  {rid}: t_end={p['t_end']}")
        if a.dry_run:
            continue
        job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=WALL,
                           template=ROOT / "config" / "slurm_mpi_template.sh", cpus=1, ntasks=NTASKS, mem="4G")
        print(f"    job={job}")


if __name__ == "__main__":
    main()
