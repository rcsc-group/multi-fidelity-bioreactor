"""Replay vs database DNS, current postprocess, read-only (task of 2026-11-04).

For each rpm of the chained replay (rc{L}_rpm{rpm}_rep) and its database DNS (kmix_l8_* for L8, fig9_l9_* for L9):
recompute kLa_1T_10 / kLa_1T_25 / kLa_1T_50 and dtmix_0.50 from the raw run files with the CURRENT
scripts/postprocess.py functions, in memory. Nothing is written into runs/ or into any results.json.
Also reports the stored results.json values (control: diary 2026-10-31 numbers came from these).

Usage: uv run python scripts/experiments/check_replay_kla25_metrics.py 8
Writes: experiments/multifidelity/check_replay_kla25_L{LEVEL}.json
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402
from scripts.experiments.replay_chain import DB, rid  # noqa: E402
from scripts.hypotheses.test_h1_mfbml import RPMS  # noqa: E402

THRESH = {"kLa_1T_10": 0.10, "kLa_1T_25": 0.25, "kLa_1T_50": 0.50}
QOIS = ("kLa_1T_10", "kLa_1T_25", "kLa_1T_50", "dtmix_0.50")


def current(run: Path) -> dict:
    params = json.loads((run / "params.json").read_text())
    t, c_star = pp._compute_c_star(run)
    T_bio, T_per = pp._t_scales(params)
    kla_to_h = 3600.0 / T_bio
    out = {q: pp._kla_period_at_threshold(t, c_star, th, T_per) * kla_to_h for q, th in THRESH.items()}
    out["dtmix_0.50"] = pp._compute_mixing_metrics(run, params)["dtmix_0.50"]
    out["span_cycles"] = float((t[-1] - t[0]) / T_per)
    out["c_star_last"] = float(c_star[-1])
    return out


def stored(run: Path) -> dict:
    f = run / "results.json"
    d = json.loads(f.read_text()) if f.exists() else {}
    return {q: d.get(q, math.nan) for q in QOIS}


def pct(a: float, b: float) -> float:
    return 100.0 * (a / b - 1.0) if (a == a and b == b and b != 0) else math.nan


def main(level: int) -> None:
    rows, control = [], []
    for rpm in RPMS:
        dns = ROOT / "runs" / DB[level].format(rpm)
        rep = rg.RUNS / rid(level, rpm, "rep")
        if not (rep / "tr_oxy.dat").exists():
            print(f"{rpm:5g} rpm  replay missing")
            continue
        d_cur, r_cur = current(dns), current(rep)
        d_sto, r_sto = stored(dns), stored(rep)
        dev = {q: pct(r_cur[q], d_cur[q]) for q in QOIS}
        rows.append({"rpm": rpm, "dev_pct": dev, "dns": d_cur, "replay": r_cur})
        control.append({"rpm": rpm,
                        "dns_recomputed_vs_stored_pct": {q: pct(d_cur[q], d_sto[q]) for q in QOIS},
                        "replay_recomputed_vs_stored_pct": {q: pct(r_cur[q], r_sto[q]) for q in QOIS},
                        "stored_dev_pct": {q: pct(r_sto[q], d_sto[q]) for q in QOIS}})
        print(f"{rpm:5g} rpm  " + "  ".join(f"{q} {dev[q]:+6.1f}" for q in QOIS)
              + f"  span {r_cur['span_cycles']:.1f} cyc")
    print("counts within 10 %:", {q: sum(1 for r in rows if abs(r["dev_pct"][q]) <= 10) for q in QOIS}, "of", len(rows))
    out = ROOT / "experiments" / "multifidelity" / f"check_replay_kla25_L{level}.json"
    out.write_text(json.dumps({"level": level, "rows": rows, "control": control}, indent=1, default=float))
    print("wrote", out)


if __name__ == "__main__":
    main(int(sys.argv[1]))
