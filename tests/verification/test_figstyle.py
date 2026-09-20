"""The replica figures share one visual grammar; these tests pin it.

The failure this guards against is not a crash -- it is two figures that each
render fine and contradict each other, which no amount of running the scripts
will reveal. #CC79A7 meant L9 in Figs 9-12 and L10 in Fig 13 for weeks. So the
assertions here are about the TABLES, plus a source-level check that no script
hard-codes a hue that the tables own.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from scripts import figstyle as fs

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"

# Figures that draw more than one mesh level, or Kim alongside us, and so
# must speak the grammar. Schematics (Fig 1) and field maps are exempt.
GRAMMAR_SCRIPTS = [
    "plot_fig2.py", "plot_fig8.py", "plot_fig9.py", "plot_fig10.py",
    "plot_fig11_fig12.py", "plot_fig13.py", "plot_fig14_fig15.py",
    "plot_figA16_c.py", "plot_figA16_current.py", "plot_figA18.py",
]


def test_kim_is_black():
    assert fs.KIM == "black"


def test_every_level_has_a_distinct_colour():
    hues = list(fs.LEVEL_COLOUR.values())
    assert len(hues) == len(set(hues))


def test_level_hues_are_disjoint_from_component_and_tracer_hues():
    """A reader who learns blue = L8 must not meet blue = u_x elsewhere."""
    levels = set(fs.LEVEL_COLOUR.values())
    assert not levels & set(fs.COMPONENT_COLOUR.values())
    assert not levels & set(fs.TRACER_COLOUR.values())
    assert fs.KIM not in levels


def test_quantity_markers_are_unique():
    """No shape may mean two different physical quantities.

    The threshold families are exempt by design -- chi = 0.50 and C* = 50%
    share "v" precisely so the loosest criterion reads the same in Fig 9 and
    Fig 11 -- and that sharing is pinned separately below.
    """
    shapes = [v for k, v in fs.MK.items()
              if not k.startswith(("chi_", "cstar_"))]
    assert len(shapes) == len(set(shapes))


def test_each_threshold_family_is_internally_unique():
    for prefix in ("chi_", "cstar_"):
        shapes = [v for k, v in fs.MK.items() if k.startswith(prefix)]
        assert len(shapes) == len(set(shapes)) == 3


def test_threshold_shapes_do_not_collide_with_quantity_shapes():
    quantities = {fs.MK[k] for k in
                  ("tau", "ediss", "vorticity", "elevation", "l2")}
    thresholds = {v for k, v in fs.MK.items()
                  if k.startswith(("chi_", "cstar_"))}
    assert not quantities & thresholds


def test_threshold_families_are_ordered_consistently():
    """Loosest criterion shares a shape across the chi and C* families, so
    the ranking reads the same in Fig 9 and Fig 11."""
    assert fs.threshold_marker("chi", 0.50) == fs.threshold_marker("cstar", 50)
    assert fs.threshold_marker("chi", 0.75) == fs.threshold_marker("cstar", 25)
    assert fs.threshold_marker("chi", 0.95) == fs.threshold_marker("cstar", 10)


def test_fill_encodes_statistic_not_provenance():
    assert fs.series_kw("k", "o", stat="max")["mfc"] == "k"
    assert fs.series_kw("k", "o", stat="mean")["mfc"] == "w"


def test_linestyle_encodes_provenance_not_statistic():
    """Kim solid, ours dotted -- and that must not change with the statistic,
    which is the overload Fig 13 used to carry."""
    assert fs.series_kw("k", "o", stat="max")["ls"] == "-"
    assert fs.series_kw("k", "o", stat="mean")["ls"] == "-"
    assert fs.series_kw("k", "o", stat="max", ours=True)["ls"] == ":"
    assert fs.series_kw("k", "o", stat="mean", ours=True)["ls"] == ":"


def test_level_colour_rejects_an_unassigned_level():
    with pytest.raises(KeyError):
        fs.level_colour(11)


def test_stat_must_be_named_explicitly():
    with pytest.raises(ValueError):
        fs.series_kw("k", "o", stat="average")


_OWNED = ({h.lower() for h in fs.LEVEL_COLOUR.values()}
          | {h.lower() for h in fs.COMPONENT_COLOUR.values()}
          | {h.lower() for h in fs.TRACER_COLOUR.values()})
_HEX = re.compile(r"#[0-9A-Fa-f]{6}")
_DOCSTRING = re.compile(r'""".*?"""|\'\'\'.*?\'\'\'', re.S)


def _code_only(src: str) -> str:
    """Source with comments and docstrings removed.

    Both are where a script EXPLAINS a hue it no longer uses -- the note in
    plot_figA16_current.py recording that its curves were crimson until the
    shared grammar arrived is exactly the kind of provenance worth keeping,
    and scanning it as if it were code would force us to delete the history
    to keep the test green.
    """
    body = _DOCSTRING.sub("", src)
    return "\n".join(ln for ln in body.splitlines()
                     if not ln.lstrip().startswith("#"))
# Named colours the old scripts used for levels or for Kim. Any reappearance
# is a script that has drifted back out of the grammar.
_BANNED_NAMES = ("royalblue", "orchid", "seagreen", "mediumorchid",
                 "crimson", "tab:orange")


@pytest.mark.parametrize("name", GRAMMAR_SCRIPTS)
def test_scripts_do_not_hardcode_owned_hues(name):
    """Hues the tables own must be referenced through figstyle, not retyped.

    A retyped hex is how the tables and the figures drift apart: the table
    says L9 is green, the script still says pink, and both are 'correct'
    in isolation.
    """
    path = SCRIPTS / name
    if not path.exists():
        pytest.skip(f"{name} not present")
    body = _code_only(path.read_text())
    offenders = {h.lower() for h in _HEX.findall(body)} & _OWNED
    assert not offenders, f"{name} hard-codes {sorted(offenders)}"


@pytest.mark.parametrize("name", GRAMMAR_SCRIPTS)
def test_scripts_do_not_use_the_old_ad_hoc_names(name):
    path = SCRIPTS / name
    if not path.exists():
        pytest.skip(f"{name} not present")
    body = _code_only(path.read_text())
    found = [c for c in _BANNED_NAMES if c in body]
    assert not found, f"{name} still uses {found}"
