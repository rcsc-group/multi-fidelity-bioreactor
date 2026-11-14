"""Fixtures of the guard tests. The guards read a finished feature bundle (tools/bundle); they never run the solver.

    MFB_BUNDLE=/path/to/work/bundle.json uv run pytest tests/bundle

Without MFB_BUNDLE every guard is skipped, with one reason. With it, a missing file, or a bundle built from
other solver sources than src/ (its src_hash differs), is an error: a guard must never pass on stale evidence.
A case that failed in the bundle fails the guards that need it, with the case's own error.
"""
import json
import os
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.bundle.bundle import cases_hash, src_hash  # noqa: E402
from tools.bundle import cases as bc  # noqa: E402

BARS = ROOT / "tools" / "bundle" / "bars.yaml"


def load_bars() -> dict:
    return yaml.safe_load(BARS.read_text())


def case(bundle: dict, name: str) -> dict:
    """The record of case `name`; fails the test with the case's recorded error if it did not run cleanly."""
    c = bundle["cases"].get(name)
    if c is None:
        pytest.fail(f"case {name} is not in the bundle")
    expected = "refused" if getattr(bc.CASES.get(name), "expect_refusal", False) else "ok"
    if c["status"] != expected:
        pytest.fail(f"case {name} {c['status']} (expected {expected}): {c['error']}")
    return c


NOT_REQUESTED = "feature bundle not requested (set MFB_BUNDLE=path/to/bundle.json; see tools/bundle/bundle.py)"


def pytest_collection_modifyitems(config, items):
    """Skip every guard with a marker, not from the fixture: pytest then prints the reason once per file,
    not once per test."""
    if not os.environ.get("MFB_BUNDLE"):
        for item in items:
            if "bundle" in getattr(item, "fixturenames", ()):
                item.add_marker(pytest.mark.skip(reason=NOT_REQUESTED))


@pytest.fixture(scope="session")
def bundle() -> dict:
    where = os.environ.get("MFB_BUNDLE")
    if not where:
        pytest.skip(NOT_REQUESTED)
    p = Path(where)
    if not p.is_file():
        pytest.fail(f"MFB_BUNDLE={p}: no such file", pytrace=False)
    b = json.loads(p.read_text())
    now = src_hash(Path(os.environ.get("MFB_BUNDLE_SRC", ROOT / "src")))   # MFB_BUNDLE_SRC: power runs (tools/bundle/mutants.py)
    if b.get("src_hash") != now:
        pytest.fail(f"stale bundle: built from src_hash {b.get('src_hash')}, src/ is now {now}. "
                    f"Run the bundle again on this source (tools/bundle/bundle.sbatch).", pytrace=False)
    if b.get("cases_hash") != cases_hash():
        pytest.fail(f"stale bundle: built from other case definitions (tools/bundle/cases.py hash "
                    f"{b.get('cases_hash')}, now {cases_hash()}). Run the bundle again.", pytrace=False)
    if b.get("fidelity") is not None:
        pytest.fail(f"bundle built with --fidelity {b['fidelity']} (a level check); the guards' bars are for "
                    f"the L6 bundle.", pytrace=False)
    mpi = p.with_name("mpi4.json")   # written by the separate MPI job (tools/bundle/bundle_mpi.sbatch)
    if mpi.is_file():
        m = json.loads(mpi.read_text())
        if m.get("src_hash") != now or m.get("cases_hash") != cases_hash():
            pytest.fail(f"stale MPI result {mpi}: built from src_hash {m.get('src_hash')}, src/ is now {now}.",
                        pytrace=False)
        b["cases"]["mpi4"] = m["case"]
    return b


@pytest.fixture(scope="session")
def bars() -> dict:
    return load_bars()
