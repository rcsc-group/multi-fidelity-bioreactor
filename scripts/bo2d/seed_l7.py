"""Submit cold L7 runs at every L8 site that has both kLa_1T_10 and tau_95 (chain arm seed, 2026-11-07)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import submit as S  # noqa: E402
C = S.C
k, t = C.load_db("kLa_1T_10"), C.load_db("tau_95")
s = k[k.level == 8].merge(t[t.level == 8][["run"]], on="run")[["rpm", "theta"]].drop_duplicates()
go = "--go" in sys.argv
print(len(s), "sites")
for r, th in s.to_numpy(float):
    S.cold(7, round(float(r), 3), round(float(th), 3), arm="seed", go=go)
