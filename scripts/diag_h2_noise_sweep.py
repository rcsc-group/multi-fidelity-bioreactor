"""H2 wiggliness: is it the assumed HF noise, or the L7 signal itself?

(1) Re-run the H2 scoring with the HF noise inflated: 5% (as run), 20%, 50%.
(2) Data analysis of L7: how much of the L7 rpm variation is signal? Compared:
    the spread of the L7 values across rpm, versus the run-to-run scatter at
    one condition (warm L7 vs cold L7 at 32.5 rpm, the only repeat available).
"""
import json, math, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import scripts.diag_mf_early_late as m              # noqa: E402
from scripts.diag_mf_l7_l9_kim import val, RPMS, TARGETS  # noqa: E402

for key in TARGETS:
    a = np.array([q for q in ((r, val(f"kmix_l7_rpm{r:g}", key), val(f"fig9_l9_rpm{r:g}", key))
                              for r in RPMS) if all(np.isfinite(q))], float)
    x, lf, hf = a.T
    line = f"{key:11s}"
    for noise in (0.05, 0.20, 0.50):
        m.REL_NOISE = noise
        res = {}
        for mid in x[1:-1]:
            for k, v in m.one_case(x, lf, hf, np.isin(x, [x[0], mid, x[-1]])).items():
                res.setdefault(k, []).append(v[0])
        line += f" | noise {noise:.0%}: MF {np.mean(res['MF']):.0%} GP {np.mean(res['HF-only GP']):.0%} ratio {np.mean(res['LF x ratio']):.0%}"
    cv7 = np.std(lf) / np.mean(lf)
    rep = abs(val("kmix_l7_rpm32.5", key) - val("fig9_ladder_l7_rpm32.5", key)) / val("kmix_l7_rpm32.5", key)
    cv9 = np.std(hf) / np.mean(hf)
    print(line)
    print(f"{'':11s}   L7 spread across rpm {cv7:.0%}, L7 repeat diff at 32.5 {rep:.0%}, L9 spread {cv9:.0%}")
