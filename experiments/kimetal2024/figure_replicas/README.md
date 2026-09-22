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
| 8(a) | τ and EDR against time, 3 cycles | `fig8_hist_l10` | **do not present — transient, see below** |
| 8(b,c) | τ and EDR distributions at the peak instant | `fig8_hist_l10` | good |
| 9 | mixing time against rpm | 9 of 10 L9 points | good shape; two outliers |
| 10 | mixing time against angle | 2 of 6 L9 points | **very partial** |
| 11 | kLa against rpm | same 9 points | **not grid-converged** |
| 12 | kLa against angle | 2 of 6 L9 points | **very partial** |
| 13 | τ and EDR against rpm and angle | L6/L8/L9/L10 | strongest agreement in the set |
| A.16(a,b) | velocity convergence with mesh | L6/L8/L10 | good |

## Caveats worth saying out loud

**Fig 9 reads much better at nine points than it did at six, and the earlier
verdict here was too harsh.** With only the high-rpm half present, L9 appeared
to cross Kim's curves; with 15–37.5 rpm in hand it tracks them in parallel
over most of the range. dt₉₅/Kim by rpm:

    15    17.5   22.5   25     27.5   30     32.5   35     37.5
    0.72  0.74   0.66   1.05   1.15   0.67   0.70   0.95   1.53

Most points sit at 0.66–0.74, a consistent offset rather than scatter. Three
do not: 25 and 27.5 rpm jump to ~1.1 while their neighbours stay near 0.7, and
37.5 reaches 1.53. Those are the open question, not the whole curve. 20 rpm is
still running.

**Fig 11 is the weak one, and it is a convergence problem, not a sampling
one.** The L6→L7→L8 column at 32.5 rpm runs 10.1× → 5.1× → 2.2× Kim, still
falling with refinement, so L9 at 1.13× is probably not converged either. Its
rpm dependence is also erratic where Kim's is smooth (2.96, 3.50, 1.96, 3.39,
3.70, 1.35, 1.13, 1.53, 1.56 against his values). Do not present as settled.

**Figs 10 and 12 are two of six angle points** (θ = 7 and 6) with θ = 5, 4, 3,
2 queued. Two points do not make a trend; they are here because the scripts
run and the numbers are reasonable (dt₅₀/Kim = 0.64 and 0.58; kLa₂₅ ratios
1.14 and 1.07), not because the figures are ready.

**Fig 8(a) was drawn on the wrong axis until 2026-09-21.** Kim's panel is a
time evolution over three unfolded cycles (t/T_p = 80–83) with both quantities
on one axes, τ left in blue and EDR right in red. Ours folded x mod 1, split
the two quantities into separate panels, and coloured by refinement level —
three departures at once. Now corrected. The cycle NUMBERS still differ and
are not forced to match: `fig8_hist_l10` continues a converged state at
t/T_p = 47, not 80, and for a settled periodic flow any three consecutive
settled cycles are equivalent; relabelling them 80–83 would assert a
correspondence with Kim's tracer protocol that this run does not have.

**Fig 8 is not settled, and that explains the rest of it.** Kim's ⟨τ'_w⟩(t) is
very nearly a pure sinusoid; ours is visibly distorted. Half-period
antisymmetry forbids EVEN harmonics in a settled signed wall shear, and
`fig8_hist_l10` has a large one that DECAYS through the run — 2×/1× = 0.880,
0.440, 0.228 over cycles 0–4, 4–7, 7–11 after its restart, halving every ~3.5
cycles and still 0.228 at the end. So the waveform is the run relaxing after
its warm start, not L10's shape, and the frames the figure uses are all
contaminated. `fig8_hist_l10b` re-records with 20 cycles (~15 to settle, the
tail to record).

Two claims made earlier and withdrawn: that the phase-fold scatter was genuine
cycle-to-cycle variability of the flow (it is this transient), and that L8
provided a clean sinusoidal control (L8 samples exactly 5.00 frames/cycle,
Nyquist 2.5, so it is structurally blind to a 2nd harmonic — its higher
harmonics fit as nonsense, 5×/1× = 501.8, which is how the rank deficiency
shows). Only a recording denser than ~8 frames/cycle can measure waveform
shape at all, and until `fig8_hist_l10b` lands we have exactly one such
recording and it is the unsettled one.

**Two quantitative gaps in Fig 8, probably the same cause.** Our τ amplitude matches (±2.2e-3
vs his ±2.0e-3 Pa) and EDR shows his two-peaks-per-cycle correctly, but our
EDR peaks near 0.68 W/m³ against his 0.40 — a factor 1.7 — and the shape is
spikier where his is smooth. Separately, our τ *distribution* at the peak
instant spans ±0.49 Pa while his histogram fits inside ±0.015: the means agree
while our tails are ~30× wider, and 7.5% of samples fall outside his window.
Both are most likely the settling transient above rather than a separate
defect, and the re-record tests that directly: if EDR falls to ~0.4 and the
tails narrow once the even harmonic is gone, there is nothing else to chase.
If they do not, a correct wall mean hiding wide tails is the grid-scale-noise
signature this project has been caught by before, and diagnosing it needs a
volume metric, not another wall number.

**Panel (a) keeps colour on the refinement level**, against Kim, who spends it
on the quantity (blue τ, red EDR) because he draws one dataset. Colour meaning
level is a project-wide convention a reader carries across all eighteen
figures, and it is worth more than matching one panel's palette; quantity is
separated by axis and line style instead — solid left is τ, dashed right is
EDR. Levels whose recordings cover different cycles simply do not appear.

**Panel (a) reads shear_stress.dat, not the frames.** Kim's panel is a time
series of spatial means, which is exactly what that file logs; the frames
exist for panels (b)/(c), which need per-cell fields, and are sparse in time
by design because each is ~10 MB. Switching source took L10 from 13.6 to 30.7
samples per cycle and L8 from 5 to a smooth sinusoid — L8's triangular look
was entirely its frame cadence, not its physics, which also settles that L8
looks smooth because it is smooth at this resolution.

The file is NON-DIMENSIONAL and needs the same ρU² / ρU³L⁻¹ factors the
frames path applies; omitting them rescaled τ by 6× and was caught only
because the peak moved from 2.2e-3 to 3e-4 Pa when the source changed.

**`diagnostic_Fig8_phasefold_*.png` are not replicas.** They fold all cycles
onto one period to show cycle-to-cycle repeatability, which Kim's panel does
not ask. Kept because they answer a question worth asking. In them: L8's
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
