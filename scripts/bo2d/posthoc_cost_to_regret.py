"""POST HOC metric (2026-11-08, chosen after held-out #1 failed the 5 % rule; user decision): core-h to first reach
regret <= R for R in 15/10/7 %, per arm and seed, from the synth_bench logs given on the command line."""
import re
import sys
from collections import defaultdict
rx = re.compile(r"^(chain|l9only)\s+seed (\d+) round \d+ core-h\s+(\d+) regret\s+([-\d.]+)% feasible True")
first = defaultdict(dict)
for f in sys.argv[1:]:
    for line in open(f):
        m = rx.match(line)
        if not m:
            continue
        arm, seed, ch, reg = m.group(1), int(m.group(2)), float(m.group(3)), float(m.group(4))
        for R in (15, 10, 7):
            if reg <= R and (arm, seed, R) not in first:
                first[(arm, seed, R)] = ch
seeds = sorted({k[1] for k in first} | {int(s) for s in ()})
for R in (15, 10, 7):
    print(f"regret <= {R}%")
    for s in sorted({k[1] for k in first}):
        c, l = first.get(("chain", s, R)), first.get(("l9only", s, R))
        ratio = f"{l / c:4.1f}x" if c and l else "  -"
        print(f"  seed {s:2d} chain {c} l9only {l} ratio {ratio}")
