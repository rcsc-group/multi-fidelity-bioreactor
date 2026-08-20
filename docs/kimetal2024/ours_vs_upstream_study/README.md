# Ours vs upstream (Kim et al.) comparison study

Investigation into the ~3-4x gap between our fork's shear-stress/EDR
results and Kim et al.'s published Fig. 8. Full narrative and dated
findings are in `diary.md` (repo root); this folder collects the
presentation-ready outputs.

## Contents

| File | What it shows |
|---|---|
| `01_star_max_shear_stress_ours_L10.mp4` | Our own fork, L10, τ field with a star marking the per-frame argmax(τ) location. No upstream equivalent exists at this time window (t=20.65-22.45). |
| `02_percentile_sensitivity_upstream.png` | Upstream's 99th/99.9th/100th percentile of \|τ\| across its 3 real settled snapshots -- shows the sharp blowup at the top of the tail. |
| `03_summary_numbers_table.png` | Pointwise |u|/τ diff by percentile, apples-to-apples data. P100's τ column (~2e3%) is the known cut-cell singularity (grows unboundedly with resolution, documented separately in `diary.md`), not a residual disagreement -- P50/P90/P99 are the numbers that matter. |
| `04_ours_vs_upstream_rampmatched_video.mp4` | The real answer: ours vs upstream, \|u\| and τ, both using their own real liquid mask, ramp forcing identical on both sides, 224 frames across the full run (t=0 to 13.3). |
| `05_ours_vs_upstream_rampmatched_heatmap.png` | Static version at one settled instant (t=12.7447), with a relative-error panel -- same apples-to-apples data as `04`. |

**Removed, superseded (2026-08-20):** the old 13-snapshot video and
single-instant heatmap used the mismatched ramp and a contaminated
`f>0.5`-only mask; both are fully replaced by `04`/`05` above, which
use the corrected data. Deleted rather than kept alongside, per project
convention (stale figures don't stay in presentation folders).

## Headline finding (updated 2026-08-20 -- resolved)

**Once the ramp forcing is made identical on both sides (it wasn't --
our fork had NO ramp at all, not just a shorter one) and both sides
use their own real `cs` liquid mask (not a reconstruction), ours and
upstream agree almost perfectly**: velocity relative error ~0.003-0.015%
and τ pointwise correlation ~0.999-1.000 at essentially every one of
224 matched snapshots spanning the full run (`t=0` to `13.3`). Raw
shear-stress sign agreement looks like only ~60-68% at first glance,
but that's fully explained: it's dominated by cells where `|tau|` is at
the numerical noise floor (~1e-8) and sign is physically meaningless
there -- restricting to the top 10% of cells by magnitude (the ones
that matter for any percentile-based shear-stress KPI) gives 93%
agreement, and the top 1% gives 100%. See `06_ours_vs_upstream_
rampmatched_video.mp4` and `diary.md`, 2026-08-20.

**The vortex thread is closed**: it was a mask-reconstruction artifact.
With the real `cs` column now available on our side too, argmax(|tau|)
in our own field moves substantially from snapshot to snapshot, as
expected for a real flow -- it does not sit still. The earlier
"stationary vortex" was very likely a fixed region the old analytic
mask reconstruction (`|y|<b_nd`, used only because our fork didn't dump
`cs` yet) misclassified as liquid.

**What this does and doesn't prove**: this validates that our fork
correctly reproduces Kim et al.'s own reference driver code, once the
ramp/mask confounds are removed -- the fork's numerics were not the
source of prior disagreement. It does not, by itself, confirm the
reference driver's output matches the published Fig. 8 numbers -- that
upstream-driver-vs-paper question is separate and not rechecked here.

See `diary.md` entries dated 2026-08-15 through 2026-08-20 for full
derivations, job IDs, and the mechanism checks that were ruled out along
the way (restart sensitivity, metric correction, sampling window, ramp
duration, resolution mismatch, AMR, liquid-mask contamination, bubble/
droplet removal).
