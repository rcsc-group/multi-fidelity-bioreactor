"""Why did kLa_1T_10 max 2sigma/mean rise (0.85 -> 1.16) after the first campaign L9 batch (2026-11-06)?

Hypothesis: one new L9 point disagrees with the L8 trend, the fitted length scale shrinks, and the band widens
between data points. Test: refit with all 4 new L9 rows dropped (expect about 0.85 if nothing else changed), then
with each one dropped alone. The point whose removal brings max 2s/mean back down is the driver.
Also lists the FAILED-GATE rows that the fit drops. Needs a compute job.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import campaign as C  # noqa: E402

QOI = sys.argv[1] if len(sys.argv) > 1 else "kLa_1T_10"
gates = C.gate_all_runs()
d_all = C.load_db(QOI)
d = C.exclude_failed_gate(d_all, gates)
dropped = d_all[~d_all["run"].isin(d["run"])]
print("FAILED-GATE rows dropped:")
print(dropped[["run", "level", "rpm", "theta", "y"]].to_string() if "y" in dropped else dropped.to_string())
new = d[d["run"].str.startswith("c2d_l9_")]
print("\nnew L9 rows:")
print(new.to_string())


def run(dd, label):
    m, info = C.fit_qoi_info(dd, QOI)
    r = C.max_crit(m)
    print(f"{label:40s} max 2s/mean {r['max']:.3f} at rpm {r['rpm']:.2f} theta {r['theta']:.2f}  n={len(dd)}")


run(d, "all data")
run(d[~d["run"].isin(new["run"])], "without all 4 new L9")
for rn in new["run"]:
    run(d[d["run"] != rn], f"without {rn}")
