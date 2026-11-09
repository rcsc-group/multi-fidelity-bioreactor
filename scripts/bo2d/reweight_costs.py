"""POST HOC: re-cost synth_bench logs with measured costs (L7 5, L8 26 incl. L7, L9 376 incl. L7+L8; l9only 371).
The picks per round are in the log; the method's choices do not depend on cost (MF rule), so only the x-axis changes.
Prints core-h to first reach regret <= R (15/10/7 %) per arm/seed and the l9only/chain ratio."""
import re
import sys
from collections import defaultdict
CC = {7: 5.0, 8: 26.0, 9: 376.0}; CL = 371.0
import os; OFF = float(os.environ.get("CHAIN_OFFSET", "0"))   # seed charge for the fair-seed chain
rx = re.compile(r"^(chain|l9only)\s+seed (\d+) round \d+ core-h\s+\d+ regret\s+([-\d.]+)% feasible (\w+) .*picks L7/L8/L9 (\d+)/(\d+)/(\d+)")
spent = defaultdict(float); first = {}
for f in sys.argv[1:]:
    for line in open(f):
        m = rx.match(line)
        if not m:
            continue
        arm, seed, reg, feas = m.group(1), int(m.group(2)), float(m.group(3)), m.group(4) == "True"
        n7, n8, n9 = (int(m.group(i)) for i in (5, 6, 7))
        if (arm, seed) not in spent and arm == "chain":
            spent[(arm, seed)] = OFF
        spent[(arm, seed)] += n7 * CC[7] + n8 * CC[8] + n9 * CC[9] if arm == "chain" else n9 * CL
        for R in (15, 10, 7, 5):
            if feas and reg <= R and (arm, seed, R) not in first:
                first[(arm, seed, R)] = spent[(arm, seed)]
for R in (15, 10, 7, 5):
    print(f"regret <= {R}%")
    for s in sorted({k[1] for k in first}):
        c, l = first.get(("chain", s, R)), first.get(("l9only", s, R))
        print(f"  seed {s:2d} chain {c} l9only {l} ratio " + (f"{l / c:.1f}x" if c and l else "-"))
