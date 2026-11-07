"""Is a BO objective interior on L8 data? kLa_1T_10 alone, and kLa_1T_10 per unit tau_95 (oxygen per shear).
If the argmax sits on the domain edge (37.5 rpm / 7 deg), BO is trivial and not a demonstration."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import campaign as C  # noqa: E402
k = C.load_db("kLa_1T_10"); t = C.load_db("tau_95")
for lev in (8, 9):
    a = k[k.level == lev][["run", "rpm", "theta", "y"]].rename(columns={"y": "kla"})
    b = t[t.level == lev][["run", "y"]].rename(columns={"y": "tau"})
    m = a.merge(b, on="run"); m["kla_per_tau"] = m.kla / m.tau
    print(f"L{lev} n={len(m)}")
    print(" top kLa:", m.sort_values("kla", ascending=False).head(4)[["rpm", "theta", "kla", "tau"]].round(3).values.tolist())
    print(" top kLa/tau:", m.sort_values("kla_per_tau", ascending=False).head(4)[["rpm", "theta", "kla", "tau", "kla_per_tau"]].round(3).values.tolist())
