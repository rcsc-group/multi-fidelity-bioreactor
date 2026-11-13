"""tools/bundle/mutants.py: planted defects that the guards must catch (no solver run here)."""
from pathlib import Path

import pytest

from tools.bundle import mutants as mu


def _src(tmp_path, text):
    d = tmp_path / "src"
    d.mkdir()
    (d / "a.h").write_text(text)
    return d


def test_patch_applies_once(tmp_path):
    out = mu.patched_copy(_src(tmp_path, "x = 1;\ny = 2;\n"), tmp_path / "m", [("a.h", "y = 2;", "y = 3;")])
    assert (out / "a.h").read_text() == "x = 1;\ny = 3;\n"


@pytest.mark.parametrize("text", ["nothing here\n", "y = 2;\ny = 2;\n"])
def test_patch_must_match_exactly_once(tmp_path, text):
    with pytest.raises(ValueError, match="a.h"):
        mu.patched_copy(_src(tmp_path, text), tmp_path / "m", [("a.h", "y = 2;", "y = 3;")])


def test_every_mutant_patch_matches_the_real_source_once(tmp_path):
    root = Path(__file__).resolve().parents[1]
    for name, m in mu.MUTANTS.items():
        mu.patched_copy(root / "src", tmp_path / name, m.patches)   # raises if a patch went stale


def test_verdict_needs_every_named_guard_to_fail():
    lines = ["FAILED tests/bundle/test_restart.py::test_s9_x[cont] - AssertionError",
             "PASSED tests/bundle/test_restart.py::test_s10_y[half-cont]",
             "ERROR tests/bundle/test_replay.py::test_s1_z - Failed"]
    assert mu.verdict(lines, ["test_s9_x[cont]", "test_s1_z"]) == {"test_s9_x[cont]": "failed", "test_s1_z": "failed"}
    v = mu.verdict(lines, ["test_s10_y"])
    assert v == {"test_s10_y": "passed"} and not mu.caught(v)
    assert mu.verdict(lines, ["test_absent"]) == {"test_absent": "not run"}
