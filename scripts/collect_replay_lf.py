"""Fill experiments/multifidelity/replay_lf_l8.csv from whatever replay runs have finished (scripts/replay_lf.py).

Columns: run_id, rpm, theta, source_run, <QoIs of the looped replay>, <same QoIs of the coupled source run, prefix hf_>,
record_gb (size of the recorded period), record_wall_s, replay_wall_s (sacct elapsed, from the jobs listed in submitted.jsonl).
Usage: uv run python scripts/collect_replay_lf.py
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import replay_gate as rg  # noqa: E402
from scripts import replay_lf as lf  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "experiments" / "multifidelity" / "replay_lf_l8.csv"


def elapsed_s(job: str | None) -> float:
    if not job:
        return float("nan")
    s = subprocess.run(["sacct", "-j", job, "-X", "-n", "-P", "-o", "ElapsedRaw"], capture_output=True, text=True).stdout.split()
    return float(s[0]) if s else float("nan")


def dir_gb(d: Path) -> float:
    return sum(f.stat().st_size for f in d.glob("flux_rec.*.bin")) / 1e9 if d.exists() else float("nan")


def main() -> None:
    sub = {}
    f = lf.REC_ROOT / "submitted.jsonl"
    if f.exists():
        for line in f.read_text().splitlines():
            j = json.loads(line)
            sub[j["rep_run"]] = j   # last submission wins
    rows = []
    for rpm, th in lf.conditions():
        lid = lf.lf_id(rpm, th)
        res = rg.run_dir(lid) / "results.json"
        if not res.exists():
            continue
        r = json.loads(res.read_text())
        src = lf.source_id(rpm, th)
        h = json.loads((lf.SRC_RUNS / src / "results.json").read_text())
        j = sub.get(lid, {})
        row = dict(run_id=lid, rpm=rpm, theta=th, source_run=src)
        row.update({q: r.get(q) for q in lf.QOIS})
        row.update({f"hf_{q}": h.get(q) for q in lf.QOIS})
        row.update(record_gb=dir_gb(lf.REC_ROOT / lid), record_wall_s=elapsed_s(j.get("job_record")),
                   replay_wall_s=elapsed_s(j.get("job_replay")))
        rows.append(row)
    if not rows:
        print("no finished replay yet")
        return
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} of {len(lf.conditions())} conditions -> {OUT}")
    for r in rows:
        d = {q: 100 * (r[q] / r[f"hf_{q}"] - 1) for q in ("kLa_1T_10", "dtmix_0.50") if r[q] and r[f"hf_{q}"]}
        print(r["run_id"], {k: round(v, 1) for k, v in d.items()})


if __name__ == "__main__":
    main()
