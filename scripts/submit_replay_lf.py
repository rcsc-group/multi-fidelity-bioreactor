"""Submit the record+replay pairs of the low-fidelity L8 estimate (see scripts/replay_lf.py).

Usage: uv run python scripts/submit_replay_lf.py --validate           # kmix_l8_rpm32.5, theta 7 only
       uv run python scripts/submit_replay_lf.py --all [--dry-run]    # every condition with a cycle-230 checkpoint
       uv run python scripts/submit_replay_lf.py --only 15:7 17.5:4   # rpm:theta
Collect the CSV when the jobs end: uv run python scripts/collect_replay_lf.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import replay_lf as lf  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--only", nargs="+", help="rpm:theta")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if a.validate:
        conds = [(32.5, 7.0)]
    elif a.only:
        conds = [(float(s.split(":")[0]), float(s.split(":")[1])) for s in a.only]
    elif a.all:
        conds = lf.conditions()
    else:
        ap.error("--validate, --all or --only")
    ids = []
    for rpm, th in conds:
        try:
            info = lf.submit_pair(rpm, th, dry_run=a.dry_run)
        except (FileNotFoundError, RuntimeError) as e:
            print(f"SKIP rpm {rpm:g} theta {th:g}: {e}")
            continue
        print(json.dumps(info))
        ids.append(info)
    if not a.dry_run:
        out = lf.REC_ROOT / "submitted.jsonl"
        with open(out, "a") as fh:
            for i in ids:
                fh.write(json.dumps(i) + "\n")


if __name__ == "__main__":
    main()
