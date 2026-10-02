"""MF test (Yi KRR-LR-GPR form), fidelity = run length, L9 rpm sweep.

LF y_l(x): dtmix_0.95 predicted from the chi(t) curve up to 0.50 (all rpm).
HF y_h(x): measured dtmix_0.95, 3 training rpm (lowest, 25, highest), rest held out.
Fitted in log space (relative errors). HF noise: NOT assumed to be a physical value;
set to 1e-3 (near interpolation), stated as such.
Baselines: LF alone; GP on the 3 HF points alone.
"""
import json, math, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import scripts.diag_mf_early_late as m                       # noqa: E402
from scripts.pilot_truncated_chi import chi_curve, RPMS      # noqa: E402

m.REL_NOISE = 1e-3
rows = []
for r in RPMS:
    t, chi = chi_curve(f"fig9_l9_rpm{r:g}")
    k50 = int(np.argmax(chi >= 0.50))
    sel = (chi >= 0.25) & (np.arange(len(t)) <= k50)
    b, a = np.polyfit(t[sel], np.log(1 - chi[sel]), 1)
    pred = (math.log(0.05) - a) / b
    meas = json.load(open(ROOT / f"runs/fig9_l9_rpm{r:g}/results.json")).get("dtmix_0.95", math.nan)
    if np.isfinite(meas):
        rows.append((r, pred, meas))
x, lf, hf = np.array(rows).T
print(f"corr(log LF, log HF) = {np.corrcoef(np.log(lf), np.log(hf))[0, 1]:+.2f}")
# one_case works in linear space with relative noise; pass logs shifted positive
L, H = np.log(lf), np.log(hf)
res = {}
for mid in x[1:-1]:
    tr = np.isin(x, [x[0], mid, x[-1]])
    out = m.one_case(x, L, H, tr)
    for k in ("MF", "HF-only GP"):
        res.setdefault(k, []).append(out[k][0])
    res.setdefault("LF alone", []).append(float(np.sqrt(np.mean(((L - H) / H)[~tr] ** 2))))
for k, v in res.items():
    print(f"{k:11s} mean held-out rel. RMSE of log(dtmix): {np.mean(v):.1%}")
print("LF alone, error in dtmix itself: median |LF/HF - 1| =", f"{np.median(np.abs(lf / hf - 1)):.0%}")
