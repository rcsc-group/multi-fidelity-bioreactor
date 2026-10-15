import sys, math, json
import numpy as np
sys.path.insert(0, "/oscar/scratch/eaguerov/kp_parity")
import predict as P
from scripts import kim_convert as kc, chi_scale
rid = sys.argv[1]
fd = P.ROOT / "runs" / rid / "fields"
files = sorted(fd.glob("snap_*.bin"))
times = np.array([kc._header_t(p) for p in files])
first = None
for k, p in enumerate(files[:40]):
    if chi_scale.read_snapshot(p, ["c2"])["c2"].max() > 0: first = k; break
ini = float(f"{times[first]:.9g}")
tl = times[-1]
dt1 = 0.1519 * 7
fin1 = ini + math.floor((tl - ini) / dt1 - 1e-9) * dt1
fin1 = float(f"{fin1:.10g}")
dt2 = round((tl - times[first]) / (len(times) - 1 - first), 9)
# truncation for KP2: first plain-chi >=0.95 index + 30 (nondim chi, grid of full window)
g = kc.kim_grid(ini, dt2, tl - dt2 * 0.5)
chi, _ = P.chi_series(files, times, g, False)
fi = int(np.argmax(chi >= 0.95)) if chi.max() >= 0.95 else None
n2 = min((fi + 31) if fi is not None else len(g), len(g))
fin2 = float(f"{ini + (n2 - 1) * dt2:.10g}")
out = dict(first=first, ini=ini, dt1=dt1, fin1=fin1, dt2=dt2, fin2=fin2, n2=n2, first_idx=fi, Tper_nd=None)
print(json.dumps(out))
# advance prediction (script seconds, RPM=32.5)
for name, dt, fin in (("kp1", dt1, fin1), ("kp2", dt2, fin2)):
    r, _ = P.run(rid, ini, dt, fin)
    print(name, "pred plain %.4f strict %.4f nfiles %d" % (r["plain"], r["strict"], r["nfiles"]))
