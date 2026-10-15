"""Test S (pre-registered, diary 2026-10-15): is the L8 scatter of dtmix / kLa made by the estimator or by the flow?

S1: the two +-0.1-rpm pairs (22.4/22.6, 29.9/30.1): max |chi_a - chi_b| and late decay rates lam (fit of ln(1 - chi)
    vs t over chi 0.5-0.95). DEFINITION-DRIVEN if max|dchi| <= 0.02 in both; FLOW-DRIVEN if lam differs >= 15% in both.
S2: E0 (test_e0_scatter.analyse, same rule) on the 21 grid runs with robust estimators:
    dtmix_rob = time at which the S1 fit line reaches chi 0.95; kLa_glob = -slope of ln(1 - C*) vs t from release to
    the first C* >= 0.5. ESTIMATOR-DRIVEN if s_rob <= 0.5 s_E0; FLOW-DRIVEN if >= 0.8 s_E0.
chi(t), C*(t) rebuilt exactly as scripts/postprocess.py (_compute_mixing_metrics, _compute_c_star).
Usage: uv run python scripts/test_s_estimator.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.postprocess import (_COL_C2_LIQ_SUM, _COL_C2_LIQ_SUM2, _COL_F_LIQ, _COL_OXY_LIQ,  # noqa: E402
                                 _joined_cols, _t_scales)
from scripts.test_e0_scatter import PAIRS, analyse, pair_run  # noqa: E402

E0 = {"dtmix_0.95": 0.141, "kLa_1T_25": 0.189}


def curves(run):
    """t [s from release], chi, C* (postprocess.py definitions)."""
    d = ROOT / "runs" / run
    params = json.loads((d / "params.json").read_text())
    t, cs, cs2, ox = _joined_cols(d, "tr_oxy.dat", [_COL_C2_LIQ_SUM, _COL_C2_LIQ_SUM2, _COL_OXY_LIQ])
    _, f = _joined_cols(d, "vol_frac_interf.dat", [_COL_F_LIQ])
    fm = f.mean()
    s2 = cs2 / fm - (cs / fm) ** 2
    i0 = int(np.where(cs > 1e-10 * fm)[0][0])
    chi = np.clip(1 - s2 / s2[i0], 0, 1)
    T_bio, _ = _t_scales(params)
    ts = (t - t[i0]) * T_bio
    return ts[i0:], chi[i0:], (ox / fm)[i0:]


def late_fit(ts, chi):
    m = (chi >= 0.5) & (chi <= 0.95)
    k = np.argmax(chi >= 0.95) if (chi >= 0.95).any() else len(chi)
    m[k + 1:] = False                                   # up to the first 0.95 crossing
    b, a = np.polyfit(ts[m], np.log(1 - chi[m]), 1)    # ln(1 - chi) = a + b t
    return -b, (np.log(0.05) - a) / b                  # lam [1/s], t at chi = 0.95 [s]


def kla_glob(ts, cstar):
    k = int(np.argmax(cstar >= 0.5)) if (cstar >= 0.5).any() else len(cstar)
    b = np.polyfit(ts[:k], np.log(np.clip(1 - cstar[:k], 1e-10, None)), 1)[0]
    return -b * 3600.0                                  # 1/h


def main():
    out = {"S1": {}, "S2": {}}
    print("S1 pairs:")
    dchi_ok, lam_ok = [], []
    for a, b in PAIRS:
        ta, ca, _ = curves(pair_run(a))
        tb, cb, _ = curves(pair_run(b))
        la, t95a = late_fit(ta, ca)
        lb, t95b = late_fit(tb, cb)
        tend = max(ta[np.argmax(ca >= 0.95)], tb[np.argmax(cb >= 0.95)])
        grid = np.linspace(0, tend, 4000)
        dchi = float(np.max(np.abs(np.interp(grid, ta, ca) - np.interp(grid, tb, cb))))
        dlam = abs(la / lb - 1)
        dchi_ok.append(dchi <= 0.02)
        lam_ok.append(dlam >= 0.15)
        out["S1"][f"{a}/{b}"] = dict(max_dchi=dchi, lam=[la, lb], dlam=dlam, t95_rob=[t95a, t95b])
        print(f"  {a}/{b} rpm: max|dchi| {dchi:.3f}; lam {la:.4f} vs {lb:.4f} 1/s (diff {100 * dlam:.0f}%); "
              f"dtmix_rob {t95a:.1f} vs {t95b:.1f} s")
    s1 = "DEFINITION-DRIVEN" if all(dchi_ok) else "FLOW-DRIVEN" if all(lam_ok) else "MIXED"
    out["S1"]["reading"] = s1
    print("  S1 reading:", s1)

    ids = json.loads((ROOT / "experiments/multifidelity/test_t5_roughness.json").read_text())["run_ids"]
    grid_r = sorted(float(r) for r in ids if abs(float(r) / 0.625 - round(float(r) / 0.625)) < 1e-9
                    and 20 <= float(r) <= 32.5)
    vals = {"dtmix_rob": [], "kLa_glob": []}
    for r in grid_r:
        ts, chi, cst = curves(ids[f"{r:g}"])
        vals["dtmix_rob"].append(late_fit(ts, chi)[1])
        vals["kLa_glob"].append(kla_glob(ts, cst))
    pv = {"dtmix_rob": [], "kLa_glob": []}
    for a, b in PAIRS:
        for name, fn in (("dtmix_rob", lambda c: late_fit(c[0], c[1])[1]), ("kLa_glob", lambda c: kla_glob(c[0], c[2]))):
            pv[name].append(abs(np.log(fn(curves(pair_run(b)))) - np.log(fn(curves(pair_run(a))))))
    print("S2 robust estimators on the 21 grid runs (E0 rule):")
    for name, ref in (("dtmix_rob", "dtmix_0.95"), ("kLa_glob", "kLa_1T_25")):
        res = analyse(np.array(grid_r), np.log(vals[name]), pv[name], 2)
        ratio = res["s"] / E0[ref]
        reading = "ESTIMATOR-DRIVEN" if ratio <= 0.5 else "FLOW-DRIVEN" if ratio >= 0.8 else "MIXED"
        res.update(ratio_to_E0=ratio, reading_S2=reading, values=vals[name])
        out["S2"][name] = res
        print(f"  {name}: s {res['s']:.3f} (E0 {ref} {E0[ref]:.3f}; ratio {ratio:.2f}) rho1 {res['rho1']:+.2f} "
              f"E0-shape {res['reading']} -> {reading}")
    (ROOT / "experiments/multifidelity/test_s_estimator.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
