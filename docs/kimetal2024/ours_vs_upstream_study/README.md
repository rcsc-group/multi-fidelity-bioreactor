# Ours vs upstream (Kim et al.) comparison study

Investigation into the ~3-4x gap between our fork's shear-stress/EDR
results and Kim et al.'s published Fig. 8. Full narrative and dated
findings are in `diary.md` (repo root); this folder collects the
presentation-ready outputs.

## Contents

| File | What it shows |
|---|---|
| `01_star_max_shear_stress_ours_L10.mp4` | Our own fork, L10, τ field with a star marking the per-frame argmax(τ) location. No upstream equivalent exists at this time window (t=20.65-22.45). |
| `02_ours_vs_upstream_13snapshot_comparison.mp4` | Ours vs upstream, \|u\| and τ side by side, animated across all 13 matching L10 snapshots (t=0 to 12.76). Domain boxed, no colorbar. |
| `03_ours_vs_upstream_single_instant_heatmap.png` | Static version of the same comparison at one instant (t=12.76), with a relative-error panel. |
| `04_percentile_sensitivity_upstream.png` | Upstream's 99th/99.9th/100th percentile of \|τ\| across its 3 real settled snapshots -- shows the sharp blowup at the top of the tail. |
| `05_summary_numbers_table.png` | All key quantitative findings in one table. |

## Headline finding

Velocity's *aggregate* behavior (mean speed, bulk mean) matches well
between codebases (~0.2-4% relative error). Shear stress does not:
pointwise correlation between the two τ fields is essentially zero
(-0.15 to +0.09) at **every one of the 13 checked snapshots** spanning
the whole run, with 14-17% of cells disagreeing on the sign of τ. Any
single "bulk mean τ error" number is unstable (ranges -22% to +1989%
across snapshots) and should not be quoted as a fixed bias -- the real
finding is the decorrelation itself, not a stable offset.

**Open thread (2026-08-19):** a vortex visible in our own τ field
appears NOT to move across snapshots where the underlying flow clearly
does -- flagged as suspicious, not yet investigated.

See `diary.md` entries dated 2026-08-15 through 2026-08-19 for full
derivations, job IDs, and the mechanism checks that were ruled out
(restart sensitivity, metric correction, sampling window, ramp
duration, resolution mismatch, AMR).
