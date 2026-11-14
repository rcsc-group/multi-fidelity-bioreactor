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
    "restart_label": ["half", "cont", "caseB_warm"],
    "chaining": ["half", "cont"],
    "cross_level_warmstart": ["l5", "xl6", "cold"],
    "flux_replay": ["rec", "replay", "loop"],
}
# kill / kill_resume are left out: kill dumps when its wall-clock limit fires, so their outputs differ run to run.
WALL_CLOCK = {"logstats.dat"}   # wall-clock timings: differ between any two runs (bars.yaml noise floor)


def carry(old_bundle: Path, new_bundle: Path, mechanism: str) -> str:
    old = json.loads(Path(old_bundle).read_text())["cases"]
    new = json.loads(Path(new_bundle).read_text())["cases"]
    compared = 0
    for case in CASES[mechanism]:
        if case not in new:
            continue   # the current bundle never ran it: nothing to compare (the guards report that)
        if case not in old or old[case].get("status") != "ok":
            return f"stale: {case} not run on the old source"
        for f, h in sorted(new[case]["files"].items()):
            if f not in WALL_CLOCK and old[case]["files"].get(f) != h:
                return f"stale: {case}/{f}"
        compared += 1
    return "current" if compared else "stale: no case of this mechanism in the new bundle"


if __name__ == "__main__":
    print(carry(Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]))
