"""tau_95 in Pa (tau_nd * rho_w * (a/T_bio)^2, postprocess.py l.904) with kLa_1T_10, for the L8 and L9 rows.
Used to pick the tau_95 cap of the BO constraint (2026-11-07)."""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2])); sys.path.insert(0, str(Path(__file__).resolve().parent))
import campaign as C  # noqa: E402
from scripts.postprocess import _t_scales  # noqa: E402
k = C.load_db("kLa_1T_10"); t = C.load_db("tau_95")
m = k.merge(t[["run", "y"]], on="run", suffixes=("_kla", "_tau"))
pa = []
for r in m.run:
    p = json.loads((C.ROOT / "runs" / r / "params.json").read_text()); Tb, _ = _t_scales(p)
    a = p.get("geometry", {}).get("a", 0.25); pa.append(1e3 * (a / Tb) ** 2)
m["tau_Pa"] = m.y_tau * pa
for lev in (8, 9):
    s = m[m.level == lev].sort_values("tau_Pa")
    print(f"L{lev}:", [(round(a, 1), round(b, 1), round(c, 1), round(d, 3)) for a, b, c, d in s[["rpm", "theta", "y_kla", "tau_Pa"]].values])
