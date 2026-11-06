"""kLa band rise after L9 rows (2026-11-06), test 2: the HF posterior over (length, noise ratio) and the
split of the variance at the max-crit point into within-node (GP) and between-node (hyperparameter mix) parts.
Hypothesis: with few HF rows, the posterior moves to short length scales / high noise, so the between-node or
within-node variance rises far from the data. Falsified if the posterior and the split do not change."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import campaign as C  # noqa: E402
Q = "kLa_1T_10"
d = C.exclude_failed_gate(C.load_db(Q), C.gate_all_runs())
new = d[d["run"].str.startswith("c2d_l9_")]["run"]
P = np.array([[35.0, 4.8]])
def show(dd, lab):
    m, info = C.fit_qoi_info(dd, Q)
    hi = dd[dd.level >= 9]
    w = m.w
    lm, em = float(w @ m.l), float(w @ m.eta)
    top = np.argsort(w)[::-1][:3]
    mu, sc, nu, ww = m.predict_dist(P)
    wm = ww @ mu
    within = float((ww[:, None] * sc**2 * (nu/(nu-2))).sum(0)[0]); between = float((ww[:, None]*(mu-wm)**2).sum(0)[0])
    print(f"{lab:28s} nHF={len(hi)} L9={int((hi.level==9).sum())} std(logY)={np.std(np.log(hi.y)):.3f} "
          f"eta_lo={info['hf_eta_floor']:.3f} E[l]={lm:.3f} E[eta]={em:.3f} "
          f"top(l,eta,w)={[(round(m.l[i],2),round(m.eta[i],3),round(w[i],2)) for i in top]} "
          f"@35/4.8 within={within:.4f} between={between:.4f} crit={2*np.sqrt(np.expm1(within+between)):.3f}")
show(d, "all")
show(d[~d.run.isin(new)], "no new L9")
for r in new: show(d[d.run != r], "no " + r)
