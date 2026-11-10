"""Go/no-go for 'search cheap, confirm at L9': on the FROZEN real truth, where is the constrained optimum of each level,
and what L9 regret does the cheap-level optimum give? Constraint tau_95 <= cap at the same level (what a cheap search sees)
and at L9 (oracle constraint). Usage: uv run python scripts/bo2d/diag_level_optima.py [cap]"""
import os
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
import synth_bench as S  # noqa: E402
from scripts.bo2d import bo  # noqa: E402

cap = float(sys.argv[1]) if len(sys.argv) > 1 else 0.014
G = bo.grid(); lc = np.log(cap)
for obj in ("kla", "mix"):
    os.environ["BO_OBJ"] = obj; S.setup("real")
    k = {l: S.F[l][0](G) for l in (7, 8, 9)}; t = {l: S.F[l][1](G) for l in (7, 8, 9)}
    opt9 = k[9][t[9] <= lc].max()
    print(f"\n{obj}: L9 optimum {G[np.argmax(np.where(t[9] <= lc, k[9], -np.inf))]}")
    for l in (7, 8):
        for cl in (l, 9):
            m = t[cl] <= lc
            x = G[np.argmax(np.where(m, k[l], -np.inf))]
            i = np.argmax(np.where(m, k[l], -np.inf))
            feas9 = t[9][i] <= lc
            print(f"  L{l} objective, tau at L{cl}: argmax {x}  L9 regret {1 - np.exp(k[9][i] - opt9):6.1%}  L9-feasible {feas9}")
    # how many grid points within 5 % regret at L9 (size of the target)
    print(f"  L9 points within 5 % regret: {np.mean((t[9] <= lc) & (k[9] >= opt9 + np.log(0.95))):.1%} of grid")
