"""Diagnosis of the uncertainty growth after adding the c2d L8 runs (lead, 2026-11-05). Read-only with respect to the
fit: nothing is written except the printed table. Compute job only.

H1 binary offset: c2d L8 (2ef5b07) vs database L8 (f1c11e0) at nearest neighbours in (rpm, theta), ln difference.
H2 group offset in EDR: ln(edr) regressed on ln(tau) with a c2d indicator; the group effect vs the floor 0.10.
H3 hyperparameters: LF selected restart and HF posterior (length, noise ratio, rho), old data vs new data.
H4 curvature at theta < 3: old model (database rows only) predicts the c2d values there; error in ln.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import campaign as c  # noqa: E402

FLOOR = {"tau_95": 0.017, "kLa_1T_10": 0.095, "edr_95": 0.0085}


def nearest_pairs(d_new):
    """H1: c2d L8 rows vs nearest database L8 row (|d theta| <= 0.3, |d rpm| <= 2)."""
    c2d = d_new[(d_new.level == 8) & d_new.run.str.startswith("c2d_")]
    db = d_new[(d_new.level == 8) & ~d_new.run.str.startswith("c2d_")]
    out = []
    for r in c2d.itertuples(index=False):
        cand = db[(np.abs(db.theta - r.theta) <= 0.3) & (np.abs(db.rpm - r.rpm) <= 2.0)]
        if cand.empty:
            continue
        k = int(np.argmin(np.hypot((cand.rpm - r.rpm) / 22.5, (cand.theta - r.theta) / 5.0).values))
        out.append(np.log(r.y) - np.log(cand.iloc[k].y))
    return np.array(out)


def hyper_summary(model, info):
    """Selected LF row and the HF posterior: weighted length, noise ratio and rho (over the grid nodes)."""
    lf_sel = next(r for r in info["lf_rows"] if r["seed"] == info["selected_seed"])
    w = np.asarray(model.w, float).ravel()
    l = np.asarray(model.l, float)
    eta = np.asarray(model.eta, float)
    rho = np.asarray(model._rho, float)
    rho_w = (w[:, None] * rho).sum(0)
    imax = int(np.argmax(w))
    return dict(lf_seed=info["selected_seed"], lf_ls=lf_sel["lengthscales"], lf_noise=round(lf_sel["noise_std"], 5),
                lf_lml=round(lf_sel["lml"], 2), hf_l_mean=float((w * l).sum()), hf_eta_mean=float((w * eta).sum()),
                hf_max_node=(round(float(l[imax]), 4), round(float(eta[imax]), 5), round(float(w[imax]), 3)),
                hf_rho_mean=np.round(rho_w, 4).tolist())


def main():
    print("=== H1: c2d (2ef5b07) minus database (f1c11e0), nearest L8 neighbours, ln units ===")
    for q in ("tau_95", "kLa_1T_10"):
        d = c.load_db(q)
        diff = nearest_pairs(d)
        if diff.size:
            print(f"{q:10s} n={diff.size:3d} mean={diff.mean():+.4f} sd={diff.std(ddof=1) if diff.size > 1 else 0:.4f} "
                  f"floor={FLOOR[q]} mean/floor={diff.mean() / FLOOR[q]:+.2f} frac>0={np.mean(diff > 0):.2f}")
        else:
            print(f"{q:10s} no pairs")

    print("\n=== H2: ln(edr) = a + b ln(tau) + g * c2d ===")
    e = c.load_db("edr_95")[["run", "y"]].rename(columns={"y": "edr"})
    t = c.load_db("tau_95")[["run", "y"]].rename(columns={"y": "tau"})
    m = e.merge(t, on="run")
    X = np.column_stack([np.ones(len(m)), np.log(m.tau.values), m.run.str.startswith("c2d_").astype(float).values])
    yv = np.log(m.edr.values)
    beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
    res = yv - X @ beta
    s2 = res @ res / (len(yv) - X.shape[1])
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    print(f"n={len(yv)} (c2d {int(X[:, 2].sum())}, backfill {int(len(yv) - X[:, 2].sum())}) b={beta[1]:.3f} "
          f"group g={beta[2]:+.4f} se={se[2]:.4f}  floor={FLOOR['edr_95'] * 0 + 0.10} (task floor 0.10)  "
          f"|g|/0.10={abs(beta[2]) / 0.10:.2f}")

    print("\n=== H3 and H4: old data (database rows only) vs new data (with c2d) ===")
    for q in ("tau_95", "kLa_1T_10", "edr_95"):
        d_new = c.load_db(q)
        d_old = d_new[~d_new.run.str.startswith("c2d_")].reset_index(drop=True)
        m_old, i_old = c.fit_qoi_info(d_old, q)
        m_new, i_new = c.fit_qoi_info(d_new, q)
        print(f"-- {q}")
        print("   old:", hyper_summary(m_old, i_old))
        print("   new:", hyper_summary(m_new, i_new))
        if q in ("tau_95", "kLa_1T_10"):
            cc = d_new[d_new.run.str.startswith("c2d_") & (d_new.theta < 3.0)]
            if len(cc):
                X = cc[["rpm", "theta"]].to_numpy(float)
                pred = np.asarray(m_old.predict_mean(X), float).ravel()
                err = pred - np.log(cc.y.values)
                print(f"   H4 theta<3 (c2d n={len(cc)}): old-model error ln(pred/y): mean={err.mean():+.3f} "
                      f"rms={np.sqrt(np.mean(err ** 2)):.3f} floor={FLOOR[q]} rms/floor={np.sqrt(np.mean(err ** 2)) / FLOOR[q]:.1f}")
                for r, er in zip(cc.itertuples(index=False), err):
                    print(f"      {r.run:28s} y={r.y:.4g} err={er:+.3f}")


if __name__ == "__main__":
    main()
