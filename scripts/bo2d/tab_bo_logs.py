"""From bo_real frozen logs: L9-only evaluations (and core-h) to reach regret levels; chain picks per level."""
import re
import sys
from collections import defaultdict
rx = re.compile(r"^(chain|l9only)\s+seed (\d+) round (\d+) core-h\s+(\d+) regret\s+([-\d.]+)% feasible (\w+) .*picks L7/L8/L9 (\d+)/(\d+)/(\d+)")
jobs = dict(l.split()[2:4][::-1] and (l.split()[3], l.split()[:3]) for l in open("scripts/bo2d/state/bo_real_frozen_jobs.txt"))
for job, (obj, seed, arm) in sorted(jobs.items(), key=lambda kv: kv[1]):
    first = {}; picks = [0, 0, 0]; last = None
    for line in open(f"logs/bo_real_{job}.out"):
        m = rx.match(line)
        if not m:
            continue
        ch, reg, feas = int(m.group(4)), float(m.group(5)), m.group(6) == "True"
        picks = [p + int(m.group(i)) for p, i in zip(picks, (7, 8, 9))]; last = ch
        for lv in (30, 20, 10, 5):
            if feas and reg <= lv and lv not in first:
                first[lv] = (ch, list(picks))
    s = " ".join(f"{lv}%:{first[lv][0]}" if lv in first else f"{lv}%:-" for lv in (30, 20, 10, 5))
    print(f"{obj} {seed} {arm:6s} final core-h {last} picks L7/L8/L9 {picks} | core-h to {s}")
