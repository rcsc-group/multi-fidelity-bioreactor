"""Grade X9 (diary 2026-11-20): L9 warm start (from L8) vs L9 cold from rest. Hydrodynamic QoIs graded at 5 %
(the 2-sd-of-period-scatter widening of the bar is not applied: stricter); scalar QoIs reported only."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import postprocess as pp  # noqa: E402

S = Path("/oscar/scratch/eaguerov/mpi_runs")
ok = True
for site in ("rpm25_th4", "rpm32.5_th7"):
    c, w = pp.main(str(S / f"ft9_l9cold_{site}")), pp.main(str(S / f"ft9_l9warm_{site}"))
    for q, g in (("tau_95_qss", 1), ("edr_95", 1), ("vel_rms_qss", 1), ("kLa_1T_10", 0), ("dtmix_0.25", 0),
                 ("dtmix_0.50", 0)):
        e = 100 * (w[q] - c[q]) / c[q]
        if g:
            ok &= abs(e) <= 5
        print(f"{site:11s} {q:12s} cold {c[q]:.4g} warm {w[q]:.4g} err {e:+.1f}%{'' if g else ' (reported)'}")
print("X9", "PASS" if ok else "FAIL")
