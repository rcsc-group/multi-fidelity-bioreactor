"""Add kLa_1T_{10,25,50} to existing results.json files, touching nothing else.

The whole-period kLa (postprocess._kla_period_at_threshold, 2026-10-01) is new.
Re-running postprocess.main on old runs would recompute every KPI with today's
code, so this script only ADDS the three new keys, computed from the run's own
oxygen series (chain-joined, as postprocess does).

Usage:  uv run python scripts/backfill_kla_1T.py [--dry-run] [glob ...]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.postprocess import _compute_c_star, _kla_period_at_threshold, _t_scales  # noqa: E402


def backfill(run: Path, dry: bool) -> str:
    res, prm = run / "results.json", run / "params.json"
    if not (res.exists() and prm.exists() and (run / "tr_oxy.dat").exists()):
        return "skip"
    params = json.loads(prm.read_text())
    try:
        t, c = _compute_c_star(run)
    except (FileNotFoundError, ValueError):
        return "skip"
    T_bio, T_per = _t_scales(params)
    to_h = 3600.0 / T_bio
    new = {f"kLa_1T_{p}": _kla_period_at_threshold(t, c, p / 100, T_per) * to_h
           for p in (10, 25, 50)}
    if all(math.isnan(v) for v in new.values()):
        return "nan"
    if not dry:
        d = json.loads(res.read_text())
        d.update(new)
        res.write_text(json.dumps(d, indent=2))
    return "ok " + " ".join(f"{v:.3g}" for v in new.values())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("globs", nargs="*", default=["*"])
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    runs = sorted({p for g in a.globs for p in (ROOT / "runs").glob(g) if p.is_dir()})
    counts = {}
    for r in runs:
        s = backfill(r, a.dry_run)
        counts[s.split()[0]] = counts.get(s.split()[0], 0) + 1
        if s.startswith("ok") and any(k in r.name for k in ("fig9_", "fig10_", "fig11_")):
            print(f"  {r.name:28s} {s}")
    print(counts)


if __name__ == "__main__":
    main()
