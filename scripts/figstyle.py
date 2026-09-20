"""One visual grammar for every Kim et al. (2024) replica figure.

Before this module each script invented its own encoding, and the inventions
contradicted each other: #CC79A7 was L9 in Figs 9-12 and L10 in Fig 13, Kim
was black in Figs 10-13 and royalblue in Fig 9, marker "s" meant EDR in Fig 13
and steady-streaming vorticity in Figs 9-10, and a hollow marker meant "this
is a spatial mean" in Fig 13 but "this is the other quantity" in Figs 9-11.
A reader who learns the legend of one figure is then actively misled by the
next one, which is worse than having no convention at all.

THE GRAMMAR -- four channels, one dimension each, never overloaded:

  colour      WHOSE DATA IT IS.  Kim is black; each of our mesh levels owns
              one hue for the life of the project.  Nothing else may use
              these hues, which is why vector components have their own
              disjoint palette below.

  linestyle   WHETHER IT IS KIM'S.  Kim solid, ours dotted.  This is the one
              distinction that must survive being printed in greyscale.

  marker      WHICH QUANTITY (or which threshold of it).  One master table,
              globally unique: no shape means two things in two figures.

  fill        WHICH STATISTIC.  Filled = an extremum (absolute max, peak
              instant).  Hollow = an average (spatial mean, time mean).

Corollaries worth stating because the old scripts broke them:

  * Levels keep their hue when a figure calls them by mesh size instead of by
    level.  A.16's n_L = 2^6 IS level 6 and is drawn in level 6's orange.
  * Thresholds of one quantity (chi = 0.50/0.75/0.95, C* = 50/25/10%) share a
    shape family ordered by stringency: v -> ^ -> P.  They are distinct from
    every physical-quantity shape, so Fig 9's chi = 0.95 marker can never be
    confused with Fig 13's tau marker.
  * Statistic lives in fill ALONE, which frees linestyle for provenance.
    Fig 13 previously spent linestyle on max-vs-mean and had nothing left to
    say "this row is Kim's".

Palette is Okabe-Ito, which stays distinguishable under the common forms of
colour blindness and in greyscale reproduces as distinct lightnesses.

Usage:
    from scripts.figstyle import KIM, level_colour, MK, series_kw, rcparams
    ax.plot(x, y, **series_kw(KIM, MK["tau"], stat="max"))
    ax.plot(x, y, **series_kw(level_colour(9), MK["tau"], stat="max",
                              ours=True))
"""
from __future__ import annotations

# --------------------------------------------------------------- colour
KIM = "black"

# Level -> hue. Fig 13 is the replica that has been checked against Kim by
# eye and approved, so its assignments (L6 orange, L8 blue, L9 green,
# L10 pink) are the ones every other figure moves to match. L7 takes the
# remaining Okabe-Ito sky blue; it previously held #009E73, which Fig 13
# had already given to L9.
LEVEL_COLOUR = {
    6:  "#E69F00",   # orange
    7:  "#56B4E9",   # sky blue
    8:  "#0072B2",   # blue
    9:  "#009E73",   # bluish green
    10: "#CC79A7",   # reddish purple
}

# Vector components in single-run figures (Figs 2, A.18). Deliberately
# disjoint from LEVEL_COLOUR: a reader who has learned that blue means L8
# must not meet a blue that means u_x three figures later.
COMPONENT_COLOUR = {"x": "#D55E00", "y": "#882255", "z": "#117733"}

# Tracer configurations in Figs 5-6. Also disjoint from the level hues.
TRACER_COLOUR = {
    "c2": "#332288",   # top half
    "c1": "#AA4499",   # left half
    "c3": "#44AA99",   # circle
    "c":  "#DDCC77",   # line
}

NEUTRAL = "0.45"        # annotations, guide lines, threshold rules
GUIDE = "0.75"          # zero lines, grid emphasis


def level_colour(level: int | float) -> str:
    """Hue for a mesh level, raising rather than silently picking a default.

    A KeyError here means a figure is plotting a level the project has not
    assigned a colour to, which must be decided once in this table and not
    improvised at the call site -- that improvisation is what produced the
    L9/L10 collision this module exists to end.
    """
    return LEVEL_COLOUR[int(level)]


# --------------------------------------------------------------- marker
# Globally unique shapes. The threshold family (v ^ P) is ordered by
# stringency so a reader can rank the three without consulting the legend,
# and shares no shape with any physical quantity.
MK = {
    # physical quantities
    "tau":       "o",    # wall shear stress
    "ediss":     "s",    # energy dissipation rate
    "vorticity": "D",    # <|xi_bar|>, steady-streaming vorticity
    "elevation": "h",    # free-surface elevation (Figs 14-15)
    "l2":        "X",    # discrete L2 error (Fig A.16c)
    # mixing time, by degree of homogeneity
    "chi_0.50":  "v",
    "chi_0.75":  "^",
    "chi_0.95":  "P",
    # kLa, by the C* threshold used to fit it; ordered to match chi, so the
    # loosest criterion is "v" in both families
    "cstar_50":  "v",
    "cstar_25":  "^",
    "cstar_10":  "P",
}


def threshold_marker(kind: str, value: float) -> str:
    """Marker for a threshold series. `kind` is "chi" or "cstar"."""
    if kind == "chi":
        return MK[f"chi_{value:.2f}"]
    if kind == "cstar":
        return MK[f"cstar_{value:g}"]
    raise KeyError(f"unknown threshold family {kind!r}")


# --------------------------------------------------------------- assembly
LW_KIM, LW_OURS = 1.3, 1.1
MS_KIM, MS_OURS = 6.0, 5.0


def series_kw(colour: str, marker: str, *, stat: str = "max",
              ours: bool = False, **extra) -> dict:
    """Keyword arguments for one plotted series, per the grammar above.

    `stat` is "max" (an extremum -> filled) or "mean" (an average -> hollow).
    `ours` selects the dotted linestyle and the slightly lighter weight that
    separate our data from Kim's without spending a colour on it.
    """
    if stat not in ("max", "mean"):
        raise ValueError(f"stat must be 'max' or 'mean', got {stat!r}")
    filled = stat == "max"
    kw = {
        "color": colour,
        "marker": marker,
        "ls": ":" if ours else "-",
        "lw": LW_OURS if ours else LW_KIM,
        "ms": MS_OURS if ours else MS_KIM,
        "mfc": colour if filled else "w",
        "mec": colour,
        "mew": 1.1,
    }
    kw.update(extra)
    return kw


def rcparams() -> dict:
    """Shared typography and axis style. Kim's figures are serif with inward
    ticks; matching that costs nothing and makes side-by-side comparison read
    as one document rather than two."""
    return {
        "mathtext.fontset": "cm",
        "font.family": "serif",
        "axes.linewidth": 1.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
    }


GRID_KW = {"which": "major", "ls": ":", "alpha": 0.4}
