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

## Headline finding (retracted 2026-08-19, see diary)

**The earlier "velocity matches well" claim below does not hold up.**
Every prior statistic here was computed with a contaminated liquid
mask (~71% dead solid cells) and with 3 of the 13 snapshots
mislabeled as "ramp-settled" when they were not (upstream's real ramp
completes at `t=11.6132`, not `~9.87`). Once corrected, only 2 of the
13 snapshots are valid ours-vs-upstream comparison points at all
(`t=11.6963`, `t=12.7596`), and between those two, velocity relative
error swings 4%->68% and τ correlation flips sign. Shear-stress sign
agreement is a coin flip (48-51%) at every snapshot checked, settled
or not. Working hypothesis: a persistent phase lag in the fluid's
oscillatory *response* (not the forcing, which is identical in both
codes) between our fork (zero ramp -- see diary, `theta_max_prev`
equals `theta_max` for this run, making the ramp mechanism a no-op)
and upstream (genuine 16.25-cycle ramp) -- not yet confirmed. See
`diary.md`, 2026-08-19 (4).

**Open thread:** a vortex visible in our own τ field appears NOT to
move across snapshots where the underlying flow clearly does --
flagged as suspicious, not yet investigated. Ruled out dump/restart as
its cause (the source run never restarted). Cannot yet rule out that
it's a rendering artifact of the analytic liquid-mask reconstruction
our fork's dumps require (we don't carry a `cs` column like upstream
does) rather than a genuine flow feature.

See `diary.md` entries dated 2026-08-15 through 2026-08-19 for full
derivations, job IDs, and the mechanism checks that were ruled out
(restart sensitivity, metric correction, sampling window, ramp
duration, resolution mismatch, AMR).
