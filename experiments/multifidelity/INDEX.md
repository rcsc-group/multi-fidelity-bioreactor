# Figure index: multi-fidelity work (state of 2026-10-08)

All paths are relative to `experiments/multifidelity/` unless they start with another folder name.
Status: **USE** = shows a result we stand behind. **CONTEXT** = background or method picture.
**SUPERSEDED** = uses the old data window or old colours; do not present.

## A. mfbml (the method of Yi et al., KRR-LR-GPR, two fidelities) on our data
Code: `scripts/mfbml_local/`, plots by `scripts/plot_mf_*.py`.

| File | What it shows | Status |
|---|---|---|
| `mf_fig13a_means.png`, `mf_fig13b_means.png` | Yi-style prediction of the Kim Fig 13 means (τ and EDR) from L7/L8 to L9 | CONTEXT |
| `mf_tau_max_rpm.png`, `mf_tau_mean_rpm_unc.png` | τ against rpm, Yi-style prediction with uncertainty | CONTEXT |
| `mf_edr_max_angle.png`, `mf_edr_mean_angle.png`, `mf_edr_mean_angle_unc.png` | EDR against rocking angle | CONTEXT |
| `mf_h1_early_to_late.png` | test H1: can early-time runs predict late-time runs | CONTEXT |
| `mf_h2_l7_to_l9.png`, `mf_h2_l8_to_l9.png` | test H2: can a coarse level predict L9 | CONTEXT |
| `../kimetal2024/figure_replicas/replicated_Fig09_mf.png`, `replicated_Fig09_extrapolated.png` | mixing time (Kim Fig 9) with the Yi-style fit and the extrapolation | CONTEXT |

## B. Diagnostics that led to gcbml (single power law, grid order, noise)
| File | What it shows | Status |
|---|---|---|
| `observed_order_dtmix.png`, `observed_order_kla.png` | measured convergence order of mixing time and kLa | CONTEXT |
| `replicate_cv.png` | run-to-run noise from replicates | CONTEXT |
| `grid_ladder_32p5.png`, `spinup_test_l8_32p5.png`, `cost_per_level.png` | the grid ladder at 32.5 rpm, spin-up test, cost per level | CONTEXT |
| `chi_decay_rate_shape.png`, `chi_tail_evidence_l9.png`, `joint_rate_chi075.png` | mixing-time tail and decay-rate studies | CONTEXT |
| `h_kernel_test_chi50.png`, `h_kernel_test_chi75.png`, `h_kernel_test_chi95.png` | first test of the h-kernel (TWY2) | CONTEXT |
| `sigma_concepts.png` | the uncertainty terms of the method | CONTEXT |
| `guide/g1_levels.png` to `guide/g5_transform.png` | figures of `docs/yi_h_guide.md` | CONTEXT |

## C. gcbml fits on the hydrodynamic data (τ_mean, τ₉₅ against rpm; h → 0 prediction)
Data: `hydro_dataset_v2.json` (last 3 cycles before release, same for every run). Fits: `scripts/gcbml_hydro_fit_v3.py`.
Colours follow `scripts/figstyle.py`. Bands: dark = epistemic 95%, light = total 95%.

| File | What it shows | Status |
|---|---|---|
| `gcbml_hydro_v3_twy2_L8_10.png` | **main result**: L8–L10, current model | **USE** |
| `gcbml_hydro_v3_twy2_L6_10.png` | L6–L10: coarse levels widen the band | **USE** |
| `gcbml_hydro_v3_twy2_L6_9.png`, `hydro_convergence_v2_L6_9.png` | L6–L9 fit; successive level differences grow (8%, 17%, 18%) | **USE** |
| `gcbml_hydro_v3_twy2_p1.0_L8_10.png` … `p3_…` | order p fixed at 1, 1.5, 2, 3: the band follows p | **USE** |
| `gcbml_hydro_v3_twy2_sc0.2_…`, `sd0.2_…`, `both0.1_…`, `both0.2_…`, `n0.001_…`, `n0.003_…`, `n0.01_…` | amplitude and noise priors changed: the band does not change | **USE** (as proof) |
| `gcbml_hydro_v3_lb_L6_10.png`, `gcbml_hydro_v3_lb_L8_10.png` | Markov (Brownian) h-kernel instead of TWY2: no improvement | CONTEXT |
| `gcbml_hydro_v3_twy2_sat_L6_10.png`, `…_sat_L8_10.png`, `…_satA3_L6_10.png`, `…_satA3_L8_10.png` | saturating error shape on the real data: band wider, no gain | CONTEXT |
| `gcbml_hydro_v2_*.png` (6 files) | consistent data, old Blues palette | SUPERSEDED by v3 |
| `gcbml_hydro_fit.png`, `gcbml_hydro_fit_L6/L7/L8.png`, `hydro_convergence.png` | first fits, old data window | SUPERSEDED |

## D. gcbml synthetic tests (package repo: `/oscar/data/dharri15/eaguerov/Github/gcbml`)
| File | What it shows | Status |
|---|---|---|
| `benchmarks/preasymptotic/results/baseline.png` | round 1: width of the L6–L10 band against the L8–L10 band, per truth | **USE** |
| `benchmarks/preasymptotic/results/round2.png` | round 2: point error of each candidate against the L8–L10 error | **USE** |
| `benchmarks/preasymptotic/results/round2.md`, `round2b.md`, `round3.md` | the criteria tables (current model, saturating, two-term, satA3) | **USE** (tables, no figure for rounds 2b and 3) |
| `docs/fig/sigma_concepts.png` | uncertainty terms, package docs | CONTEXT |

## E. Kim et al. replicas (all other figures)
`../kimetal2024/figure_replicas/replicated_Fig*.png` (16 files) and its `README.md`.

## Numbers and reasoning behind C and D
`diary.md` (entries 2026-10-05 to 2026-10-07).
