# Kim et al. (2024) figure replicas

Every PNG here is regenerated from run data by a script in `scripts/`; none is
hand-edited. Figures carry no status text by design — what is provisional and
why lives in this file.

Regenerate everything currently backed by data:

```
uv run python scripts/plot_fig1.py
uv run python scripts/plot_fig2.py                      # defaults to fig9_l9_rpm32.5
uv run python scripts/plot_fig7_figA17.py fig9_l9_rpm32.5
uv run python scripts/plot_fig8.py
uv run python scripts/plot_fig9.py
uv run python scripts/plot_fig11_fig12.py
uv run python scripts/plot_fig13.py
uv run python scripts/plot_figA16_current.py
```

## Reading any of them

One visual grammar across the whole set (`scripts/figstyle.py`), so a legend
learned once holds everywhere:

| channel | means | values |
|---|---|---|
| colour | whose data | Kim **black**; L6 orange, L7 sky, L8 blue, L9 green, L10 pink |
| line style | provenance | Kim solid, ours dotted |
| marker | quantity, or which threshold of it | τ `o`, EDR `s`, vorticity `D`; thresholds `v`→`^`→`P` by stringency |
| fill | statistic | filled = an extremum, hollow = a mean |

L10 is Kim's own mesh. Levels are refinement levels: n_L = 2^L cells across
the bag, so A.16's `n_L = 2^8` is L8 and carries L8's colour.

## What is here

| figure | shows | backing data | confidence |
|---|---|---|---|
| 1 | the two frames of reference | schematic, no run | exact — every number is a parameter |
| 2 | rocking angle and mean velocities | `fig9_l9_rpm32.5` | good |
| 7, A.17 | oxygen transfer and the kLa fit | `fig9_l9_rpm32.5` | good |
| 8(a) | ⟨τ⟩ against rocking phase | L8, L9, `fig8_hist_l10` | good — see the note below |
| 8(b,c) | τ and EDR distributions at the peak instant | `fig8_hist_l10` | good |
| 9 | mixing time against rpm | 6 of 10 L9 points | **partial, and see the caveat** |
| 11 | kLa against rpm | same 6 points | **partial, not grid-converged** |
| 13 | τ and EDR against rpm and angle | L6/L8/L9/L10 | strongest agreement in the set |
| A.16(a,b) | velocity convergence with mesh | L6/L8/L10 | good |

## Caveats worth saying out loud

**Fig 9 and 11 are six of ten points and do not behave.** L9 crosses Kim's
curves rather than paralleling them: 27.5 rpm gives dt₉₅/Kim = 1.15 while
30 rpm, its neighbour, gives 0.67. Non-monotone neighbours that far apart is a
signal, not scatter. Fig 11 is worse — the L6→L7→L8 column at 32.5 rpm runs
10.1× → 5.1× → 2.2× Kim, still falling with refinement, so L9 at 1.1× is
probably not converged either. Do not present these as settled.

**Fig 8(a): L8 looks smoother than L10 and contains less information.** L8's
frames are spaced 0.2006 T_p, so they pile into five narrow arcs — largest gap
in phase coverage 0.106 against 0.007 for uniform sampling — and within an arc
consecutive points are consecutive cycles, so cycle-to-cycle variation reads as
gentle drift. L10 on the off-period cadence covers phase uniformly (gap 0.016
vs 0.013 uniform), so neighbouring points come from different cycles and the
same ~15% variability appears as scatter. The band is that spread; it is a
property of the flow, not of the mesh. L9 is still the old five-phase
recording and is drawn as five points joined by straight lines.

**Fig 8 previously reported the wrong peak.** The old L10 source sampled
exactly five phases however long it ran, and reported ⟨τ⟩ = 0.0016 Pa; the
re-recording (75 distinct phases) peaks near 0.0027. The old number was a
sample between peaks.

**A.16 does not use Kim's own curve colours.** His legend has n_L = 2¹⁰ in
grey. Matching him in this one figure meant contradicting every other, since
n_L = 2⁶ *is* level 6 and was drawn crimson while level 6 is orange elsewhere.
Axis ranges still follow his PDF exactly.

## Not here yet

| figure | needs | state |
|---|---|---|
| 3, 4, 5, 6 | `submit_figset.py`, ~28 h | not submitted |
| 10, 12 | angle sweep, 6 jobs | **running** (submitted 2026-09-21) |
| 14, 15 | `submit_fig14.py`, ~23 h | not submitted |
| A.16(c) | `submit_figA16_c.py`, ~32 h | parked by request |
| A.18 | `submit_figA18.py`, ~6 h | not submitted |

Fig 12 is deliberately absent rather than present-and-empty: its angle sweep is
still running, and an empty set of axes in a folder of results is a trap.
