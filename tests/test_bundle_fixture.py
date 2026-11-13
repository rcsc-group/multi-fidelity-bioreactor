"""tests/bundle/conftest.py: the guards read evidence, so missing or stale evidence must fail, never pass.

Each test runs a small pytest session (pytester) whose conftest re-exports the real bundle fixtures.
"""
import json
from pathlib import Path

import pytest

from tools.bundle.bundle import src_hash

pytest_plugins = ["pytester"]
ROOT = Path(__file__).resolve().parents[1]

GUARDS = """
from tests.bundle.conftest import case

def test_a(bundle):
    assert bundle["cases"]

def test_b(bundle):
    assert bundle["cases"]

def test_c(bundle):
    case(bundle, "half")
"""


@pytest.fixture
def guards(pytester):
    pytester.makeini("[pytest]\n")   # own config: a pyproject.toml in a parent of the temp dir must not apply
    pytester.makeconftest(f"import sys\nsys.path.insert(0, {str(ROOT)!r})\nfrom tests.bundle.conftest import *\n")
    pytester.makepyfile(test_guards=GUARDS)
    return pytester


def _bundle(tmp_path, h, half_status="ok", half_error=""):
    cases = {"cold": {"status": "ok", "error": ""}, "half": {"status": half_status, "error": half_error}}
    p = tmp_path / "bundle.json"
    p.write_text(json.dumps({"src_hash": h, "git_commit": None, "builds": {}, "cases": cases}))
    return p


def test_not_requested_skips_once(guards, monkeypatch):
    monkeypatch.delenv("MFB_BUNDLE", raising=False)
    r = guards.runpytest("-rs")
    r.assert_outcomes(skipped=3)
    reasons = [ln for ln in r.outlines if ln.startswith("SKIPPED [")]
    assert len(reasons) == 1, reasons          # one reason, said once
    assert "MFB_BUNDLE" in reasons[0]


def test_missing_bundle_fails(guards, monkeypatch, tmp_path):
    monkeypatch.setenv("MFB_BUNDLE", str(tmp_path / "nowhere" / "bundle.json"))
    r = guards.runpytest()
    r.assert_outcomes(errors=3)
    assert "nowhere" in r.stdout.str()


def test_stale_bundle_fails(guards, monkeypatch, tmp_path):
    monkeypatch.setenv("MFB_BUNDLE", str(_bundle(tmp_path, "0" * 64)))
    r = guards.runpytest()
    r.assert_outcomes(errors=3)
    out = r.stdout.str()
    assert "0" * 64 in out and src_hash(ROOT / "src") in out   # both hashes named


def test_failed_case_fails_with_its_error(guards, monkeypatch, tmp_path):
    p = _bundle(tmp_path, src_hash(ROOT / "src"), "failed", "exit 134: dump write error")
    monkeypatch.setenv("MFB_BUNDLE", str(p))
    r = guards.runpytest()
    r.assert_outcomes(passed=2, failed=1)
    assert "exit 134: dump write error" in r.stdout.str()


MPI_GUARD = """
def test_mpi(bundle):
    assert bundle["cases"]["mpi4"]["status"] == "ok"
"""


def _mpi_sidecar(tmp_path, h):
    (tmp_path / "mpi4.json").write_text(json.dumps({"src_hash": h, "case": {"status": "ok", "error": ""}}))


def test_mpi_sidecar_is_merged(pytester, monkeypatch, tmp_path):
    pytester.makeini("[pytest]\n")
    pytester.makeconftest(f"import sys\nsys.path.insert(0, {str(ROOT)!r})\nfrom tests.bundle.conftest import *\n")
    pytester.makepyfile(test_mpi=MPI_GUARD)
    h = src_hash(ROOT / "src")
    monkeypatch.setenv("MFB_BUNDLE", str(_bundle(tmp_path, h)))
    _mpi_sidecar(tmp_path, h)
    pytester.runpytest().assert_outcomes(passed=1)
    _mpi_sidecar(tmp_path, "1" * 64)            # MPI result from other sources: stale, an error
    r = pytester.runpytest()
    r.assert_outcomes(errors=1)
    assert "1" * 64 in r.stdout.str()


def test_bars_cover_every_case():
    from tests.bundle.conftest import load_bars
    from tools.bundle.cases import CASES
    bars = load_bars()
    assert set(bars) == set(CASES) - {"cold"}   # cold is the reference, every other case has a bar entry
    for name, b in bars.items():
        assert b.get("source"), f"bar {name} cites no source"
