"""Mark stuck 'new' tasks whose first run dir already exists as 'dup' (2026-11-10 21:55; see loop.propose).
Prints each task it changes. Dry run unless --go.
  uv run python scripts/bo2d/mark_dup.py [--go]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import loop as L  # noqa: E402

go = "--go" in sys.argv
tasks = L.load_tasks()
for t in tasks:
    if t["state"] != "new" or t["i"] != 0:
        continue
    step = t["steps"][0]
    lev = 8 if step == "spin8" else int(step[1])
    rid = L.C.run_id(lev, *t["x"], tag=L.TAGS[t["arm"]])
    if (L.C.ROOT / "runs" / rid).exists():
        print("dup", t["arm"], t["obj"], t["x"], step, rid)
        t["state"] = "dup"
if go:
    L.save_tasks(tasks)
    print("saved")
