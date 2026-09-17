"""Feed the Fig 9 L9 sweep into the 64-CPU per-user cap, one batch at a time.

These jobs run on the `priority` QOS, whose MaxTRESPU is cpu=312 (NOT the
`normal` QOS's cpu=64 -- always read the limit for the QOS the jobs actually
use: `squeue -o %q`). Ten 32-rank points is 320 CPUs, 8 over the cap, and on
2026-09-16 the scheduler CANCELLED nine of them outright ("CANCELLED by 0")
rather than queueing them. 9 concurrent 32-rank jobs (288 CPUs) is the most
that fits.

This picks the next rpm that has neither a results.json nor a live job, and
submits it only if doing so keeps this user's bioreactor CPU total at or under
the cap. Idempotent: safe to run on a timer.

It also finishes points that ran out of clock. t_end is sized from Kim's own
dtmix_0.95 and the first completed point overshot it (chi topped out at 0.90),
so a point can complete cleanly and still report NaN for the series Fig 9 is
about. scripts/autoextend.py continues such a point from its own checkpoint;
extensions are submitted BEFORE new points because they are short and they
close a hole rather than opening one.

Usage:  uv run python scripts/fig9_sweep_driver.py [--cap 64] [--dry-run]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
RANKS = 32

# One entry per sweep this driver feeds. `values` are submitted highest-first
# (the cheapest points at the top of the range finish soonest and free the
# queue); `target_key` is the KPI a point must reach before it counts as done,
# which differs between the figures -- Kim publishes chi=0.95 for the rpm
# series and only chi=0.50 for the angle series.
SWEEPS = {
    "fig9": {
        "values": [37.5, 35, 32.5, 30, 27.5, 25, 22.5, 20, 17.5, 15],
        "run_id": "fig9_l9_rpm{v:g}",
        "target_key": "dtmix_0.95", "target_chi": 0.95,
        "submit": ["scripts/submit_fig9.py", "--stage", "sweep", "--rpms"],
        "rpm_of": lambda v: v,
    },
    "fig10": {
        "values": [7.0, 6.0, 5.0, 4.0, 3.0, 2.0],
        "run_id": "fig10_l9_th{v:g}",
        "target_key": "dtmix_0.50", "target_chi": 0.50,
        "submit": ["scripts/submit_fig10.py", "--angles"],
        "rpm_of": lambda v: 32.5,
    },
}


def _squeue_cpus() -> tuple[int, set[str]]:
    """(CPUs held by this user's bioreactor jobs, run names in flight)."""
    out = subprocess.run(
        ["squeue", "-u", "eaguerov", "-h", "-o", "%C|%j|%T"],
        capture_output=True, text=True).stdout
    cpus, names = 0, set()
    for line in out.splitlines():
        if not line.strip():
            continue
        c, name, _state = line.split("|")
        if "bioreactor" in name:
            cpus += int(c)
            names.add(name)
    return cpus, names


def _live_run_ids() -> set[str]:
    """run_ids with a job in the queue.

    Matched via each job's PARAMS env var (the sbatch --export carries the
    params.json path, which contains the run_id). An earlier version grepped
    `scontrol show job` for the run name and always returned empty, because
    WorkDir is the repo root -- the run name never appears there. That silently
    made every in-flight point look un-submitted.
    """
    live = set()
    jids = subprocess.run(["squeue", "-u", "eaguerov", "-h", "-o", "%i"],
                          capture_output=True, text=True).stdout.split()
    for jid in jids:
        log = ROOT / "logs" / f"slurm_{jid}.out"
        text = log.read_text(errors="replace") if log.exists() else ""
        if not text:
            # PENDING jobs have no log yet -- fall back to the submission record
            rec = ROOT / "logs" / f"submitted_{jid}.txt"
            text = rec.read_text(errors="replace") if rec.exists() else ""
        for sweep in SWEEPS.values():
            for v in sweep["values"]:
                name = sweep["run_id"].format(v=v)
                if name in text:
                    live.add(name)
    return live


def _extend_short_points(sweep: dict, cpus: int, cap: int,
                         live: set[str], dry: bool) -> int:
    """Continue any point that stopped short of the target. Returns CPUs spent."""
    from scripts.autoextend import plan_extension, submit_extension, walltime_safety
    from scripts.cost_model import min_per_cycle

    spent = 0
    for v in sweep["values"]:
        base = sweep["run_id"].format(v=v)
        if base in live:
            continue
        plan = plan_extension(ROOT / "runs", base,
                              target_key=sweep["target_key"],
                              target_chi=sweep["target_chi"])
        if plan is None:
            continue
        if plan.get("blocked"):
            print(f"  {base}: BLOCKED -- {plan['reason']}")
            continue
        if plan["run_id"] in live:
            continue
        if cpus + spent + RANKS > cap:
            print(f"  {plan['run_id']}: no headroom for an extension")
            continue
        mpc, conf = min_per_cycle(level=9, ntasks=RANKS, rpm=sweep["rpm_of"](v))
        submit_extension(plan, ntasks=RANKS, min_per_cycle=mpc,
                         walltime_safety=walltime_safety(conf), dry=dry)
        if not dry:
            spent += RANKS
    return spent


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=312)  # priority QOS MaxTRESPU
    ap.add_argument("--sweep", choices=sorted(SWEEPS), default="fig9")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    sweep = SWEEPS[a.sweep]
    cpus, _ = _squeue_cpus()
    live = _live_run_ids()
    cpus += _extend_short_points(sweep, cpus, a.cap, live, a.dry_run)
    todo = [v for v in sweep["values"]
            if not (ROOT / "runs" / sweep["run_id"].format(v=v)
                    / "results.json").exists()
            and sweep["run_id"].format(v=v) not in live]

    print(f"bioreactor CPUs in flight: {cpus}/{a.cap} | live: {sorted(live)}")
    print(f"still to run: {todo}")
    if not todo:
        print("sweep complete or fully in flight")
        return
    if cpus + RANKS > a.cap:
        print(f"no headroom ({cpus} + {RANKS} > {a.cap}); nothing submitted")
        return

    v = todo[0]
    cmd = [sys.executable, *sweep["submit"], f"{v:g}"]
    if a.dry_run:
        cmd.append("--dry-run")
    print(f"submitting {v:g} ...")
    subprocess.run(cmd, cwd=ROOT)


if __name__ == "__main__":
    main()
