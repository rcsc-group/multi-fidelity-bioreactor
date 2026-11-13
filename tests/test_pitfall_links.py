"""Every `pitfall "<anchor>"` written in the code must name a section of docs_site/developer/solver-pitfalls.md.

The page is where the silent failure modes are written down; the pointer in the code is how a reader finds
it at the moment they touch that code. A pointer to a missing section is a lesson lost.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs_site" / "developer" / "solver-pitfalls.md"
POINTER = re.compile(r'[Pp]itfall "([a-z0-9-]+)"')


def _anchors() -> set:
    return {m.group(1) for m in re.finditer(r"^## ([a-z0-9-]+)\s*$", PAGE.read_text(), re.M)}


def _pointers() -> dict:
    out = {}
    for pattern in ("src/*.c", "src/*.h", "scripts/**/*.py", "tools/**/*.py"):
        for f in ROOT.glob(pattern):
            for m in POINTER.finditer(f.read_text(errors="ignore")):
                out.setdefault(m.group(1), []).append(str(f.relative_to(ROOT)))
    return out


def test_every_pitfall_pointer_has_a_section():
    missing = {a: files for a, files in _pointers().items() if a not in _anchors()}
    assert not missing, f"pitfall pointers with no section in {PAGE.name}: {missing}"


def test_the_restart_lessons_are_pointed_to_from_the_code():
    """Found by the feature bundle 2026-11-13; each must be findable from the code that can repeat it."""
    p = _pointers()
    for anchor in ("restart-ramp", "restart-phase", "record-landing"):
        assert anchor in p, f"no code points to pitfall {anchor}"
