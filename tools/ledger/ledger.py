"""Feature ledger: what each solver/workflow mechanism saves (core-h) and what error it costs.

Numbers come from the production pairs in pairs.yaml (clean only when provenance.pair_status says
same-source) and from the L6 feature bundle (tools/bundle). Cost is read from sacct: every job whose submit
line names the run's params.json, summed (chained segments, timeouts and failures included).
"""
import math
import re
import json
import subprocess
from pathlib import Path

DTMIX_CACHE = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/experiments/multifidelity/dtmix_cache.json")
_PARAMS = re.compile(r"PARAMS=\S*/([^/,\s]+)/params\.json")


def parse_sacct(rows: list[str]) -> dict:
    """{run name: [CPU seconds of each job]} from `sacct -P -o JobID,WorkDir,NCPUS,ElapsedRaw,CPUTimeRAW,State,
    SubmitLine` rows. Jobs whose submit line names no params.json are skipped."""
    out: dict[str, list[int]] = {}
    for row in rows:
        f = row.split("|")
        if len(f) < 7:
            continue
        m = _PARAMS.search("|".join(f[6:]))
        if m:
            out.setdefault(m.group(1), []).append(int(f[4] or 0))
    return out


def sacct_jobs(since: str = "2026-08-01") -> dict:
    """parse_sacct over all of this user's jobs since `since` (run it in a batch job, not on a login node)."""
    r = subprocess.run(["sacct", "-X", "-n", "-P", "-S", since, "-o",
                        "JobID,WorkDir,NCPUS,ElapsedRaw,CPUTimeRAW,State,SubmitLine%2000"],
                       capture_output=True, text=True, check=True)
    return parse_sacct(r.stdout.splitlines())


def cost_core_h(run: str, jobs: dict) -> float | None:
    """Core-hours of every job of run `run` (None: no job found, e.g. the run predates sacct retention)."""
    s = jobs.get(run)
    return sum(s) / 3600 if s else None


def _qoi(run_dir: Path, q: str) -> float:
    f = Path(run_dir) / "results.json"
    v = json.loads(f.read_text()).get(q) if f.is_file() else None
    if v is None and q.startswith("dtmix_") and DTMIX_CACHE.is_file():   # older runs: dtmix only in the cache
        v = (json.loads(DTMIX_CACHE.read_text()).get(Path(run_dir).name) or {}).get(q)
    return float(v) if isinstance(v, (int, float)) and v is not None else math.nan


def rel_error(a: Path, b: Path, qoi: str) -> float:
    """a / b - 1 for QoI `qoi` of two run folders (NaN if either is missing or NaN)."""
    x, y = _qoi(a, qoi), _qoi(b, qoi)
    return x / y - 1 if math.isfinite(x) and math.isfinite(y) and y != 0 else math.nan


def _median(xs: list) -> float | None:
    xs = sorted(x for x in xs if x is not None and math.isfinite(x))
    if not xs:
        return None
    n = len(xs)
    return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def build(pairs: dict, status_of, jobs: dict, carry_of=None) -> dict:
    """One row per mechanism. status_of(with_dir, without_dir, flag) -> "same-source"|"confounded".
    error[qoi] = median signed relative error (with / without - 1) over clean pairs;
    saves_core_h = median (cost without - cost with) over clean pairs with both costs known.
    carry_of(mechanism) -> "current"|"stale: ..." (carry.py: old-source bundle vs HEAD bundle) sets the status
    of a row with clean pairs; without it such a row is "current"."""
    rows = {}
    for mech, m in pairs.items():
        runs, clean = [], []
        for p in m["pairs"]:
            st = status_of(p["with"], p["without"], m["flag"])
            runs.append({"with": p["with"].name, "without": p["without"].name, "pair": st,
                         "cost_with": cost_core_h(p["with"].name, jobs),
                         "cost_without": cost_core_h(p["without"].name, jobs),
                         "error": {q: rel_error(p["with"], p["without"], q) for q in m["qoi"]}})
            if st == "same-source":
                clean.append(runs[-1])
        err = {q: _median([r["error"][q] for r in clean]) for q in m["qoi"]} if clean else {}
        saves = _median([r["cost_without"] - r["cost_with"] for r in clean
                         if r["cost_with"] is not None and r["cost_without"] is not None])
        status = "confounded" if not clean else (carry_of(mech) if carry_of else "current")
        rows[mech] = {"status": status, "n_clean": len(clean), "error": err,
                      "saves_core_h": saves, "runs": runs}
    return rows


if __name__ == "__main__":   # python -m tools.ledger.ledger [OUT.json]   (in a batch job: it calls sacct)
    import sys
    from tools.ledger import provenance as pv
    known = pv.known_builds()
    cache: dict = {}

    def prov(d):
        if d not in cache:
            cache[d] = pv.provenance(d, known)
        return cache[d]

    import argparse
    from tools.ledger.carry import carry
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", type=Path, help="bundle.json built from the production source (e.g. 2ef5b07)")
    ap.add_argument("--new", type=Path, help="bundle.json built from HEAD's source")
    ap.add_argument("out", type=Path, nargs="?",
                    default=Path(__file__).resolve().parents[2] / "experiments" / "feature_ledger.json")
    a = ap.parse_args()

    def carry_of(mech):
        if not (a.old and a.new and a.old.is_file() and a.new.is_file()):
            return "stale: no old bundle"
        return carry(a.old, a.new, mech)

    rows = build(pv.load_pairs(), lambda w, wo, flag: pv.pair_status(prov(w), prov(wo), flag), sacct_jobs(),
                 carry_of=carry_of)
    out = a.out
    out.write_text(json.dumps(rows, indent=1, default=str) + "\n")
    for k, r in rows.items():
        print(f"{k:22s} {r['status'][:40]:40s} clean {r['n_clean']}/{len(r['runs'])}  saves {r['saves_core_h']}  "
              + "  ".join(f"{q} {v:+.3f}" if v is not None else f"{q} -" for q, v in r["error"].items()))
