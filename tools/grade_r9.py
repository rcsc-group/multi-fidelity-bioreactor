"""Grade R9 (diary 2026-11-20): replay vs coupled record run, kLa_1T_10 and dtmix_0.50 within 5 % at both sites."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import postprocess as pp  # noqa: E402

S = Path("/oscar/scratch/eaguerov/mpi_runs")
ok = True
SITES = sys.argv[1:] or ["rpm32.5_th7", "rpm37.5_th7"]   # e.g. rpm37.5_th7_k8
for site in SITES:
    a, b = pp.main(str(S / f"ft9_rec_{site}")), pp.main(str(S / f"ft9_rep_{site}"))
    for q in ("kLa_1T_10", "dtmix_0.50", "dtmix_0.95"):
        e = 100 * (b[q] - a[q]) / a[q]
        gate = q in ("kLa_1T_10", "dtmix_0.50")
        ok &= (abs(e) <= 5) if gate else True
        print(f"{site} {q:11s} coupled {a[q]:.5g} replay {b[q]:.5g} err {e:+.1f}%{'' if gate else ' (reported)'}")
print("R9", "PASS" if ok else "FAIL")
