"""kLa band rise, test 3: split the within-node variance at 35 rpm / 4.8 deg into its parts, weight-averaged
over the HF posterior: sigma2 = b_n/a_n * ystd^2 (fitted residual scale), v_gp (1 - k'K^-1 k), v_trend
(linear-basis uncertainty), v_noise (eta^2). Which part grows when the L9 rows are added?"""
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
    m, _ = C.fit_qoi_info(dd, Q)
    Xn = m.normalize_input(P); ks = m._kern(m._d2(m.xs, Xn)); fs = m._basis(P); Kk = m._Kinv @ ks
    R = fs.T[None] - np.einsum("np,gnm->gpm", m.F, Kk)
    vgp = (1.0 - np.einsum("gnm,gnm->gm", ks, Kk))[:, 0]
    vtr = np.einsum("gpm,gpq,gqm->gm", R, m._Ainv, R)[:, 0]
    s2 = m.b_n / m.a_n * m.yh_std**2 * m.nu / (m.nu - 2)
    w = m.w
    f = lambda a: float(w @ a)
    print(f"{lab:26s} ystd={m.yh_std:.3f} sigma2={f(s2):.3f} vgp={f(vgp):.3f} vtrend={f(vtr):.3f} veta={f(m.eta**2):.3f} "
          f"s2*vgp={f(s2*vgp):.3f} s2*vtr={f(s2*vtr):.3f} s2*eta2={f(s2*m.eta**2):.3f} nu={float(np.ravel(m.nu)[0]):.1f}")
    print("   HF rows:", dd[dd.level >= 9][["rpm", "theta", "y"]].round(2).values.tolist())
show(d, "all"); show(d[~d.run.isin(new)], "no new L9")
