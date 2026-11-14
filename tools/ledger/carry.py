"""Does a mechanism's evidence, measured on an old solver source, still hold for the current source?

    carry.py OLD_BUNDLE NEW_BUNDLE MECHANISM   ->   "current" | "stale: <first differing output>"

Both arguments are bundle.json files (tools/bundle) built from the two sources. The mechanism's cases (CASES
below) must have produced byte-identical outputs on both (the bundle stores each file's sha256). Bitwise is
stricter than the guard bars: a harmless round-off change marks the row stale, and the old pair is then
re-measured, never silently trusted.
"""
import json
import sys
from pathlib import Path

CASES = {   # bundle cases whose outputs carry each mechanism (tools/bundle/cases.py)
    "restart_label": ["half", "cont", "kill", "kill_resume", "caseB_warm"],
    "chaining": ["half", "cont", "kill", "kill_resume"],
    "cross_level_warmstart": ["l5", "xl6", "cold"],
    "flux_replay": ["rec", "replay", "loop"],
}


def carry(old_bundle: Path, new_bundle: Path, mechanism: str) -> str:
    old = json.loads(Path(old_bundle).read_text())["cases"]
    new = json.loads(Path(new_bundle).read_text())["cases"]
    for case in CASES[mechanism]:
        if case not in new:
            continue   # the current bundle never ran it: nothing to compare (the guards report that)
        if case not in old or old[case].get("status") != "ok":
            return f"stale: {case} not run on the old source"
        for f, h in sorted(new[case]["files"].items()):
            if old[case]["files"].get(f) != h:
                return f"stale: {case}/{f}"
    return "current"


if __name__ == "__main__":
    print(carry(Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]))
