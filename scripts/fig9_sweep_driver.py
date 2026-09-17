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
RPMS = [37.5, 35, 32.5, 30, 27.5, 25, 22.5, 20, 17.5, 15]
RANKS = 32


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
        for rpm in RPMS:
            name = f"fig9_l9_rpm{rpm:g}"
            if name in text:
                live.add(name)
    return live


def _extend_short_points(cpus: int, cap: int, live: set[str], dry: bool) -> int:
    """Continue any point that stopped short of chi=0.95. Returns CPUs spent."""
    from scripts.autoextend import plan_extension, submit_extension, walltime_safety
    from scripts.cost_model import min_per_cycle

    spent = 0
    for rpm in RPMS:
        base = f"fig9_l9_rpm{rpm:g}"
        if base in live:
            continue
        plan = plan_extension(ROOT / "runs", base)
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
        mpc, conf = min_per_cycle(level=9, ntasks=RANKS, rpm=rpm)
        submit_extension(plan, ntasks=RANKS, min_per_cycle=mpc,
                         walltime_safety=walltime_safety(conf), dry=dry)
        if not dry:
            spent += RANKS
    return spent


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=312)  # priority QOS MaxTRESPU
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    cpus, _ = _squeue_cpus()
    live = _live_run_ids()
    cpus += _extend_short_points(cpus, a.cap, live, a.dry_run)
    todo = [r for r in RPMS
            if not (ROOT / "runs" / f"fig9_l9_rpm{r:g}" / "results.json").exists()
            and f"fig9_l9_rpm{r:g}" not in live]

    print(f"bioreactor CPUs in flight: {cpus}/{a.cap} | live: {sorted(live)}")
    print(f"still to run: {todo}")
    if not todo:
        print("sweep complete or fully in flight")
        return
    if cpus + RANKS > a.cap:
        print(f"no headroom ({cpus} + {RANKS} > {a.cap}); nothing submitted")
        return

    rpm = todo[0]
    cmd = [sys.executable, "scripts/submit_fig9.py", "--stage", "sweep",
           "--rpms", f"{rpm:g}"]
    if a.dry_run:
        cmd.append("--dry-run")
    print(f"submitting {rpm:g} rpm ...")
    subprocess.run(cmd, cwd=ROOT)


if __name__ == "__main__":
    main()
