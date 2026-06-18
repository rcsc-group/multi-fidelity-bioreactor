"""Check pass/fail status of MPI checkpoint-restart grid runs.

Usage:
    uv run python scripts/check_restart_run.py runs/grid_H*
    uv run python scripts/check_restart_run.py runs/grid_H12 runs/grid_H14

Exit 0 if all specified runs pass, 1 otherwise.

Pass criteria:
  1. results.json exists with finite kLa_25 and kLa_50
  2. tr_oxy.dat exists with at least 1 row of finite oxy_liq_sum > 0
  3. No 'SIGFPE' or 'SIGSEGV' or 'Killed' in SLURM error log
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path


def check_run(run_dir: Path) -> tuple[bool, str]:
    """Returns (pass, detail_string)."""
    issues = []

    # Check results.json
    results_file = run_dir / "results.json"
    if results_file.exists():
        try:
            r = json.loads(results_file.read_text())
            kla25 = r.get("kLa_25", None)
            kla50 = r.get("kLa_50", None)
            if kla25 is None or not math.isfinite(float(kla25)):
                issues.append(f"kLa_25 not finite: {kla25}")
            if kla50 is None or not math.isfinite(float(kla50)):
                issues.append(f"kLa_50 not finite: {kla50}")
        except Exception as e:
            issues.append(f"results.json parse error: {e}")
    else:
        issues.append("results.json missing")

    # Check tr_oxy.dat (canonical run dir first, then scratch fallback)
    scratch_run = Path("/oscar/scratch/eaguerov/mpi_runs") / run_dir.name
    tr_oxy_canonical = run_dir / "tr_oxy.dat"
    tr_oxy_scratch = scratch_run / "tr_oxy.dat"
    oxy_path = tr_oxy_canonical if tr_oxy_canonical.exists() else (
        tr_oxy_scratch if tr_oxy_scratch.exists() else None
    )

    if oxy_path:
        lines = [l for l in oxy_path.read_text().splitlines() if l.strip()]
        if not lines:
            issues.append("tr_oxy.dat empty")
        else:
            last = lines[-1].split()
            try:
                oxy_val = float(last[2])
                if not math.isfinite(oxy_val) or oxy_val <= 0:
                    issues.append(f"oxy_liq_sum not physical: {oxy_val}")
            except (IndexError, ValueError) as e:
                issues.append(f"tr_oxy.dat parse error: {e}")
    else:
        issues.append("tr_oxy.dat not found (may still be running)")

    # Check SLURM error logs for crash signals
    crash_keywords = ["SIGFPE", "SIGSEGV", "Killed", "Segmentation fault", "Floating point exception"]
    err_files = []
    for log_root in [run_dir / "logs", scratch_run / "logs"]:
        if log_root.exists():
            err_files += list(log_root.glob("slurm_*.err"))
    for ef in err_files:
        content = ef.read_text(errors="replace")
        for kw in crash_keywords:
            if kw in content:
                issues.append(f"crash signal '{kw}' in {ef.name}")
                break

    passed = len(issues) == 0
    detail = "; ".join(issues) if issues else "PASS"
    return passed, detail


def main():
    if len(sys.argv) < 2:
        print("Usage: check_restart_run.py runs/grid_H*", file=sys.stderr)
        sys.exit(1)

    run_dirs = [Path(a) for a in sys.argv[1:]]
    all_passed = True
    results = []

    for rd in run_dirs:
        if not rd.exists():
            print(f"  {rd.name:25s}  NOT FOUND")
            all_passed = False
            results.append((rd.name, False, "directory not found"))
            continue
        passed, detail = check_run(rd)
        status = "PASS" if passed else "FAIL"
        print(f"  {rd.name:25s}  {status}  {detail}")
        if not passed:
            all_passed = False
        results.append((rd.name, passed, detail))

    print()
    n_pass = sum(1 for _, p, _ in results if p)
    print(f"{n_pass}/{len(results)} runs passed")
    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
