"""Static checks on the solver source (fast, no build).

A bare `return;` inside a Basilisk event returns an undefined int, and a nonzero value stops the time loop
(docs_site/developer/solver-pitfalls.md#bare-return). This test fails if any event body in src/ has one.
"""
import re
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
EVENT = re.compile(r"^\s*event\s+\w+\s*\(.*\)\s*\{?\s*$")


def event_bodies(text: str):
    """Yield (event header line number, body text) for every event with a body."""
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        if EVENT.match(lines[i]) and not lines[i].rstrip().endswith(";"):
            start, depth, body = i, 0, []
            j = i
            opened = False
            while j < len(lines):
                code = lines[j].split("//")[0]
                depth += code.count("{") - code.count("}")
                opened = opened or "{" in code
                body.append(code)
                if opened and depth == 0:
                    break
                j += 1
            yield start + 1, "\n".join(body)
            i = j
        i += 1


def test_finds_events():
    """The parser must see the events (guards against a test that silently checks nothing)."""
    n = sum(1 for f in SRC.glob("*.[ch]") for _ in event_bodies(f.read_text()))
    assert n >= 15, n


def test_no_bare_return_in_events():
    bad = [f"{f.name}:{ln}" for f in sorted(SRC.glob("*.[ch]"))
           for ln, body in event_bodies(f.read_text()) if re.search(r"\breturn\s*;", body)]
    assert not bad, f"bare `return;` in event(s): {bad}"
