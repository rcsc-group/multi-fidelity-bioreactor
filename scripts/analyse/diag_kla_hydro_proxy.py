"""Can a hydrodynamics-only scalar stand in for kLa as Yi et al.'s LF output?

Hypothesis: kLa = kL * a, and both factors are set by the settled flow alone
(oxygen is passive). Two classical closures need no scalar transport:
  eddy-cell (Lamont & Scott):  kL ~ Sc^-1/2 (eps nu)^1/4  ->  P_eddy = a eps^1/4
  penetration (Higbie):        kL ~ sqrt(D u / l)         ->  P_pen  = a sqrt(u)
with a = interface length (interface_area(f), 2D) and eps, u taken from the
QSS window BEFORE release, i.e. from what a hydrodynamics-only run would have.

Falsifier: if log P does not track log kLa25 across the 14 L9 points (rank
correlation well below ~0.7), the proxy is not an informative LF and Yi's
linear transfer has nothing to transfer.
"""
import glob, json, math, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr

R = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs/"
KIM = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/experiments/kimetal2024/csv_raw/"
kim_f = pd.read_csv(KIM + "mixing_kla_vs_frequency.csv").set_index("RPM")
kim_a = pd.read_csv(KIM + "mixing_kla_vs_angle.csv").set_index("Angle_deg")
kcol = "kLa_exp5pts_25"

rows = []
runs = [r for r in sorted(glob.glob(R + "fig9_l9_rpm*")) + sorted(glob.glob(R + "fig10_l9_th*"))
        if not r.endswith("_ext1") and os.path.exists(r + "/results.json")]
for d in runs:
    p = json.load(open(d + "/params.json")); res = json.load(open(d + "/results.json"))
    rpm = round(p["omega_b"] * 60 / (2 * math.pi), 2); th = p["theta_max"][0]
    T = 2 * math.pi / p["omega_b"]; L, H = .25, .0715
    V = L / 4 * (H + .5 * L * math.tan(math.radians(th))); U = V / (H * .5) / T
    Tnd = T / (L / U)
    t_inj = (p.get("t_checkpoint") or 0) + p["n_mix_cycles"] * Tnd
    lo = t_inj - 20 * Tnd
    vf = pd.read_csv(d + "/vol_frac_interf.dat", sep=r"\s+")
    m = (vf.t > lo) & (vf.t < t_inj)
    a = vf.f_liq_interf[m].mean() * L                        # m (per unit depth)
    eps = res["ediss_mean_qss"]                              # W/m3
    u = res["vel_rms_qss"] * U                               # m/s
    kim = kim_f.loc[rpm, kcol] if d.find("fig9") >= 0 else kim_a.loc[th, kcol]
    rows.append(dict(run=os.path.basename(d), rpm=rpm, th=th, a=a, eps=eps, u=u,
                     kla=res["kLa_25"], kim=kim, n=int(m.sum())))
df = pd.DataFrame(rows)
df["P_eddy"] = df.a * df.eps ** 0.25
df["P_pen"] = df.a * np.sqrt(df.u)
pd.set_option("display.width", 200)
print(df.round(4).to_string(index=False))
for P in ("P_eddy", "P_pen", "a", "eps", "u"):
    for tgt in ("kla", "kim"):
        rs = spearmanr(df[P], df[tgt])[0]; rp = pearsonr(np.log(df[P]), np.log(df[tgt]))[0]
        print(f"{P:7s} vs {tgt:4s}: spearman {rs:+.2f}  pearson(log) {rp:+.2f}")
