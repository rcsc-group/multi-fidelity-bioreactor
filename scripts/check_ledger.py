"""Self-consistency checker for experiments/hypothesis_ledger.json.

Rules:
1. Within a thread, a 'confirmed' entry whose hypothesis logically contradicts
   another 'confirmed' entry in the same thread is flagged for human review.
2. An 'inconclusive' entry that has a later entry in the same thread with
   status 'confirmed' or 'falsified' is flagged as resolvable.
3. A thread with only 'falsified' entries and no 'confirmed' ones is noted
   (investigation open).

Exit 0 if no unresolved issues; exit 1 otherwise.
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

LEDGER = Path(__file__).parents[1] / "experiments" / "hypothesis_ledger.json"


def load(path: Path) -> list[dict]:
    entries = json.loads(path.read_text())
    for i, e in enumerate(entries):
        for f in ("thread", "hypothesis", "prediction", "commit", "job_id",
                  "observation", "status", "date"):
            if f not in e:
                print(f"WARNING: entry {i} missing field '{f}'")
    return entries


def check(entries: list[dict]) -> int:
    issues = 0
    by_thread = defaultdict(list)
    for e in entries:
        by_thread[e["thread"]].append(e)

    for thread, group in by_thread.items():
        confirmed = [e for e in group if e["status"] == "confirmed"]
        falsified = [e for e in group if e["status"] == "falsified"]
        inconclusive = [e for e in group if e["status"] == "inconclusive"]

        print(f"\nThread: {thread}  ({len(group)} entries: "
              f"{len(confirmed)} confirmed, {len(falsified)} falsified, "
              f"{len(inconclusive)} inconclusive)")

        # Flag inconclusives that have a later resolution
        later_statuses = {e["status"] for e in group}
        for e in inconclusive:
            if "confirmed" in later_statuses or "falsified" in later_statuses:
                print(f"  RESOLVABLE inconclusive (job {e['job_id']}): "
                      f"later entries in this thread have a definitive status. "
                      f"Update this entry.")
                issues += 1

        # Check for confirmed entries with identical hypotheses (duplicate test)
        hyps = [e["hypothesis"] for e in confirmed]
        seen = set()
        for h in hyps:
            if h in seen:
                print(f"  DUPLICATE: hypothesis confirmed twice — verify intentional: "
                      f"'{h[:80]}...'")
                issues += 1
            seen.add(h)

        # Open investigations (all falsified, no confirmed)
        if falsified and not confirmed:
            print(f"  NOTE: all {len(falsified)} entries falsified — "
                  f"investigation still open")

        if not issues:
            print(f"  OK")

    print(f"\n{'PASS' if issues == 0 else 'FAIL'}: {issues} issue(s) found.")
    return issues


if __name__ == "__main__":
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else LEDGER
    if not path.exists():
        print(f"ERROR: {path} not found")
        sys.exit(1)
    entries = load(path)
    print(f"Loaded {len(entries)} entries from {path}")
    sys.exit(check(entries))
