"""Production L9 warm-start runs: slow-flow (streaming vorticity) drift after the tracer release (diary 2026-11-21).
For each run: streaming at the release period, at the end, and the last-5-period slope. A run whose streaming still
moves after release mixes in a flow that is not its settled one (X9: -42 % dtmix_0.50 at 22 % streaming deficit)."""
import json
import sys
from pathlib import Path
import numpy as np

S = Path("/oscar/scratch/eaguerov/mpi_runs")
runs = sorted(S.glob("c2d_l9_*")) + sorted(S.glob("ft9_l9warm_*"))
names = {r.name for r in runs}
print(f"{'run':30s} {'n_rel':>5s} {'n_per':>5s} {'s_rel':>6s} {'s_end':>6s} {'drift%':>7s} {'slope%/per':>10s}")
for r in runs:
    if r.name + "_fix" in names:
        continue                                     # the _fix rerun supersedes it
    f = r / "streaming.dat"
    if not f.exists():
        continue
    s = np.atleast_2d(np.loadtxt(f, skiprows=1))[:, 2]
    n = json.loads((r / "params.json").read_text()).get("n_mix_cycles", 0)
    rel = s[min(n, len(s) - 1)]
    sl = 100 * np.polyfit(np.arange(5), s[-5:], 1)[0] / s[-1]
    print(f"{r.name:30s} {n:5d} {len(s):5d} {rel:6.2f} {s[-1]:6.2f} {100 * (s[-1] - rel) / s[-1]:+7.1f} {sl:+10.2f}")
