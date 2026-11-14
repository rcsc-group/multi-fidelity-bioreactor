"""tools/ledger/ledger.py and carry.py: cost of a run, error of a pair, and carry-over of a mechanism's
evidence from an old source to the current one (synthetic inputs; no sacct, no solver)."""
import json
import math
from pathlib import Path

import pytest

from tools.ledger import carry, ledger

ROWS = [   # JobID|WorkDir|NCPUS|ElapsedRaw|CPUTimeRAW|State|SubmitLine, as `sacct -P` prints them
    "100|/w|32|3600|115200|COMPLETED|sbatch --ntasks=32 --export=NONE,PARAMS=/s/mpi_runs/run_a/params.json,DUMP=x a.sh",
    "101|/w|32|1800|57600|TIMEOUT|sbatch --ntasks=32 --export=NONE,PARAMS=/s/mpi_runs/run_a/params.json a.sh",
    "102|/w|32|60|1920|FAILED|sbatch --ntasks=32 --export=NONE,PARAMS=/s/mpi_runs/run_a/params.json a.sh",
    "103|/w|4|360|1440|COMPLETED|sbatch --export=NONE,PARAMS=/s/mpi_runs/run_ab/params.json a.sh",
]


def test_cost_sums_every_segment_of_a_run():
    jobs = ledger.parse_sacct(ROWS)
    assert ledger.cost_core_h("run_a", jobs) == pytest.approx((115200 + 57600 + 1920) / 3600)
    assert ledger.cost_core_h("run_ab", jobs) == pytest.approx(0.4)   # no prefix confusion with run_a
    assert ledger.cost_core_h("never_ran", jobs) is None


def _res(d: Path, q: dict) -> Path:
    d.mkdir(parents=True)
    (d / "results.json").write_text(json.dumps(q))
    return d


def test_rel_error(tmp_path):
    a = _res(tmp_path / "a", {"dtmix_0.25": 12.0, "edr_95": float("nan")})
    b = _res(tmp_path / "b", {"dtmix_0.25": 10.0, "edr_95": 0.2})
    assert ledger.rel_error(a, b, "dtmix_0.25") == pytest.approx(0.2)
    assert math.isnan(ledger.rel_error(a, b, "edr_95"))
    assert math.isnan(ledger.rel_error(a, b, "kLa_1T_10"))   # missing QoI


def _bundle(tmp: Path, name: str, files: dict) -> Path:
    p = tmp / name / "bundle.json"
    p.parent.mkdir()
    p.write_text(json.dumps({"src_hash": name, "cases": {c: {"status": "ok", "files": f} for c, f in files.items()}}))
    return p


def test_carry_equal_outputs_is_current(tmp_path):
    f = {"xl6": {"tr_oxy.dat": "h1", "logstats.dat": "h2"}, "cold": {"tr_oxy.dat": "h3"}}
    old, new = _bundle(tmp_path, "old", f), _bundle(tmp_path, "new", f)
    assert carry.carry(old, new, "cross_level_warmstart") == "current"


def test_carry_one_changed_output_is_stale_and_named(tmp_path):
    f_old = {"xl6": {"tr_oxy.dat": "h1", "logstats.dat": "h2"}, "cold": {"tr_oxy.dat": "h3"}}
    f_new = {"xl6": {"tr_oxy.dat": "CHANGED", "logstats.dat": "h2"}, "cold": {"tr_oxy.dat": "h3"}}
    old, new = _bundle(tmp_path, "old", f_old), _bundle(tmp_path, "new", f_new)
    assert carry.carry(old, new, "cross_level_warmstart") == "stale: xl6/tr_oxy.dat"


def test_carry_case_missing_in_old_bundle_is_stale(tmp_path):
    old = _bundle(tmp_path, "old", {"cold": {"tr_oxy.dat": "h3"}})
    new = _bundle(tmp_path, "new", {"xl6": {"tr_oxy.dat": "h1"}, "cold": {"tr_oxy.dat": "h3"}})
    assert carry.carry(old, new, "cross_level_warmstart") == "stale: xl6 not run on the old source"


def test_build_uses_only_clean_pairs(tmp_path):
    """Error and saving are medians over same-source pairs; a mechanism with none is "confounded" and has no
    number, but its runs are still listed."""
    a, b = _res(tmp_path / "w1", {"q": 11.0}), _res(tmp_path / "c1", {"q": 10.0})
    c, d = _res(tmp_path / "w2", {"q": 30.0}), _res(tmp_path / "c2", {"q": 10.0})
    pairs = {"good": {"flag": None, "qoi": ["q"], "pairs": [{"with": a, "without": b}, {"with": c, "without": d}]},
             "bad": {"flag": None, "qoi": ["q"], "pairs": [{"with": c, "without": d}]}}
    status = {("w1", "c1"): "same-source", ("w2", "c2"): "confounded"}
    jobs = {"w1": [3600 * 10], "c1": [3600 * 50], "w2": [3600], "c2": [3600]}
    rows = ledger.build(pairs, lambda w, wo, flag: status[(w.name, wo.name)], jobs)
    assert rows["good"]["status"] == "current"
    assert rows["good"]["error"]["q"] == pytest.approx(0.1)
    assert rows["good"]["saves_core_h"] == pytest.approx(40.0)
    assert rows["good"]["n_clean"] == 1 and len(rows["good"]["runs"]) == 2
    assert rows["bad"]["status"] == "confounded" and rows["bad"]["error"] == {} and rows["bad"]["saves_core_h"] is None


def test_rel_error_reads_the_dtmix_cache(tmp_path, monkeypatch):
    """Older runs keep dtmix_* in the main repo's dtmix_cache.json, not in results.json."""
    a, b = _res(tmp_path / "a", {}), _res(tmp_path / "b", {"dtmix_0.25": 10.0})
    cache = tmp_path / "dtmix_cache.json"
    cache.write_text(json.dumps({"a": {"dtmix_0.25": 15.0}}))
    monkeypatch.setattr(ledger, "DTMIX_CACHE", cache)
    assert ledger.rel_error(a, b, "dtmix_0.25") == pytest.approx(0.5)


def test_build_marks_a_clean_row_stale_when_its_evidence_does_not_carry(tmp_path):
    """A clean production number measured on an old source holds for HEAD only if the mechanism's bundle outputs
    are equal on both sources (carry.py); otherwise the row is "stale: <first difference>" and keeps its number.
    A confounded row stays confounded whatever the carry says."""
    a, b = _res(tmp_path / "w1", {"q": 11.0}), _res(tmp_path / "c1", {"q": 10.0})
    pairs = {"m1": {"flag": None, "qoi": ["q"], "pairs": [{"with": a, "without": b}]},
             "m2": {"flag": None, "qoi": ["q"], "pairs": [{"with": a, "without": b}]},
             "m3": {"flag": None, "qoi": ["q"], "pairs": [{"with": a, "without": b}]}}
    status = {"m1": "same-source", "m2": "same-source", "m3": "confounded"}
    carried = {"m1": "current", "m2": "stale: cont/normf.dat", "m3": "stale: x"}
    rows = {}
    for m in pairs:
        rows |= ledger.build({m: pairs[m]}, lambda w, wo, flag, m=m: status[m], {}, carry_of=carried.get)
    assert rows["m1"]["status"] == "current"
    assert rows["m2"]["status"] == "stale: cont/normf.dat" and rows["m2"]["error"]["q"] == pytest.approx(0.1)
    assert rows["m3"]["status"] == "confounded"
