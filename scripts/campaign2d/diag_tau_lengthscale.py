"""Diagnosis of the tau_95 LF rpm length scale (lead, 2026-11-05). No model is changed; nothing is written but the table.

1. Closed-form leave-one-out residuals of the tau_95 LF fit (L8 rows, new data), in ln units, against rpm, split at
   theta < 3 and theta >= 4 (the database has no L8 below theta 4 before the c2d runs).
2. Drop-one refit: each theta < 3 L8 point removed in turn, the LF refitted with the same selection rule (LF_SEEDS, LML,
   noise floor). Reports the selected rpm length scale: which point(s) drive the short length scale.
Usage: uv run python scripts/campaign2d/diag_tau_lengthscale.py  (compute job)
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import campaign as c  # noqa: E402

QOI = "tau_95"


def lf_rows(d):
    lo = d[d.level == 8]
    X = lo[["rpm", "theta"]].to_numpy(float)
    Y = np.log(lo.y.to_numpy(float)).reshape(-1, 1)
    lf, rows, _ = c.fit_lf_selected(X, Y, c.noise_floor(QOI))
    return lf, rows, lo


def loo_residuals(lf):
    """Closed form: r_i = [K^-1 y]_i / [K^-1]_ii (scaled ln units), times y_std to ln."""
    K = lf._training_kernel_matrix(scaled_x=lf.sample_x_scaled, scaled_noise_std=lf.noise_std / lf.y_std)
    Kinv = np.linalg.inv(K)
    y = np.asarray(lf.sample_y_scaled, float).ravel()
    return (Kinv @ y) / np.diag(Kinv) * lf.y_std


def main():
    d = c.load_db(QOI)
    lf, rows, lo = lf_rows(d)
    ls_full = np.asarray(lf.kernel.param, float)
    print(f"full data: L8 n={len(lo)}, selected seed {next(r['seed'] for r in rows if r['ok'] and r['lml'] == max(x['lml'] for x in rows if x['ok']))}, "
          f"lengthscales (rpm, theta) = {np.round(ls_full, 4).tolist()}, noise {float(np.asarray(lf.noise_std).ravel()[0]):.4f}")
    res = loo_residuals(lf)
    lo = lo.reset_index(drop=True).assign(loo_ln=res)
    for name, m in (("theta < 3", lo.theta < 3.0), ("theta >= 4", lo.theta >= 4.0)):
        s = lo[m].sort_values("rpm")
        print(f"\nLOO residuals ({name}), n={len(s)}, rms={np.sqrt(np.mean(s.loo_ln ** 2)):.3f}")
        for r in s.itertuples(index=False):
            print(f"   rpm {r.rpm:6.2f}  theta {r.theta:4.2f}  run {r.run:28s} loo_ln {r.loo_ln:+.3f}")

    print("\nDrop-one refit (theta < 3 points removed one at a time):")
    cand = lo[lo.theta < 3.0]
    for r in cand.itertuples(index=False):
        keep = lo.run != r.run
        dd = lo[keep]
        X = dd[["rpm", "theta"]].to_numpy(float)
        Y = np.log(dd.y.to_numpy(float)).reshape(-1, 1)
        lf2, rows2, _ = c.fit_lf_selected(X, Y, c.noise_floor(QOI))
        print(f"   drop {r.run:28s} (rpm {r.rpm:6.2f}, theta {r.theta:4.2f}) -> lengthscales "
              f"{np.round(np.asarray(lf2.kernel.param, float), 3).tolist()}")


if __name__ == "__main__":
    main()
