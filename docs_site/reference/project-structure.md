# Project structure

```
multi-fidelity-bioreactor/
├── src/                      # Basilisk solver; map, event order, switches: src/README.md
│   ├── BioReactor.c          # start here: switches, globals, main() (setup), event order
│   ├── init.h scalars.h forcing.h snapshots.h checkpoint.h diagnostics.h   # one header per concern
│   ├── restart.h replay.h replay_events.h movies.h
│   ├── henry_oxy2.h          # two-phase scalars with Henry's law (oxygen, tracers)
│   ├── utils2.h              # statsf2(): statistics over the liquid
│   ├── params_read.h         # JSON → BioreactorParams struct (jsmn-based)
│   └── jsmn.h                # Minimal JSON parser (MIT, Serge Zaitsev)
│
├── scripts/                  # map: scripts/README.md. Top level = shared modules that other code imports
│   ├── simulate.py           # Core API: run_local() / submit_slurm() / wait_for_result()
│   ├── postprocess.py        # KPIs from run output files → results.json (see Reference → Output files)
│   ├── launch.py             # set up run directory + SLURM script (no submit)
│   ├── chain.py              # Workflow B: chained sweep (one param, YAML config)
│   ├── sweep.py              # Workflow E: JSON multi-param sweep (zip / cartesian; self-submitting chains)
│   ├── sample.py             # Workflow C: LHS / random / grid / Sobol batch sampling
│   ├── loop.py               # Workflow D: multi-fidelity BO loop
│   ├── suggest.py            # EI acquisition: suggest next HF point
│   ├── train_surrogate.py    # Train KRR-LR-GPR multi-fidelity surrogate
│   ├── stage_segment.py      # Safely stage a chain-restart checkpoint (refuses unless the
│   │                         # predecessor segment has a genuine results.json)
│   ├── health_report.py      # Numerical-health checks of a finished run
│   ├── figstyle.py           # One visual grammar for every figure
│   ├── run/                  # Scripts that submit or watch SLURM jobs (submit_*, *_watchdog)
│   ├── analyse/              # Analyses that read run data: collect_results.py, check_restart_run.py, ...
│   ├── figures/              # Figures and videos: plot_heatmaps.py, render_videos.py (called by SLURM), ...
│   ├── hypotheses/           # test_*.py: one numerical experiment per hypothesis in the lab diary
│   │                         # (NOT pytest tests; many submit jobs)
│   ├── campaign2d/, bo2d/    # The 2D optimisation campaign and its Bayesian optimisation
│   └── mfbml_local/          # Multi-fidelity surrogate (KRR-LR-GPR) and Bayesian variants
│
├── config/
│   ├── slurm_template.sh, slurm_mpi_template.sh, slurm_video_template.sh, slurm_health_template.sh
│   │                              # SBATCH templates for each binary variant
│   ├── param_space.yaml          # Parameter bounds for optimization suite
│   ├── bo_config.yaml            # Workflow D config (multi-fidelity BO)
│   ├── sample_config.yaml        # Workflow C config (batch sampling)
│   ├── chain_config.yaml         # Workflow B config (generic sweep template)
│   ├── chain_config_smoke.yaml   # Workflow B config (fidelity-3 smoke test)
│   ├── sweep_example.json        # Workflow E example / smoke test (2-segment, fidelity 3)
│   └── sweep_fb_theta*.json, sweep_fb_fill*.json  # theta/fill × omega_b sweeps at various fidelities
│
├── tests/
│   ├── conftest.py               # Shared fixtures (run_bioreactor, CANONICAL_PARAMS); autouse
│   │                              # session fixture rebuilds build/build-video before tests run
│   ├── test_chain.py             # chain.py unit tests (build_chain)
│   ├── test_sweep.py             # sweep.py unit tests (detect, expand, group, build, self-submit)
│   ├── test_sample.py            # sample.py unit tests (row_to_params, run_sampling)
│   ├── test_param_schema.py      # params_read schema tests
│   ├── test_launch.py            # launch.py unit tests
│   ├── test_simulate.py          # simulate.py unit tests (mocked sbatch)
│   ├── test_suggest.py           # suggest.py / EI acquisition unit tests
│   ├── test_train_surrogate.py   # surrogate training unit tests
│   ├── test_loop.py               # loop.py unit tests (wall-clock time tracking, ExperimentData)
│   ├── test_loop_integration.py   # end-to-end MF-BO testbed on a synthetic objective (see
│   │                               # Workflow D → Testing the algorithm cheaply)
│   ├── test_stage_segment.py     # stage_segment.py checkpoint-staging safety tests
│   ├── test_postprocess.py       # postprocessing unit tests
│   ├── test_shear_stress.py      # shear-stress percentile KPI pipeline tests
│   ├── test_health_report.py     # health_report.py numerical-health-check tests
│   ├── integration/
│   │   └── test_bioreactor_runs.py  # end-to-end output file checks (real binary, `medium`)
│   └── verification/
│       ├── test_grid_convergence.py   # velocity RMS convergence L5 vs L6 (`hpc`)
│       ├── test_mass_conservation.py, test_oxygen_monotonicity.py, test_forcing_frequency.py
│       └── ...                        # other numerical-verification tests (all `medium`)
│
├── examples/
│   └── synthetic_bo_demo.py  # Runnable companion to "Your first optimization loop" -- no
│                              # Basilisk build or SLURM access needed
├── experiments/
│   └── figures/              # Heatmap PDFs generated by plot_heatmaps.py (committed)
├── docs/                     # Internal project notes (specs, literature, meetings, solver history)
│                              # -- NOT this documentation site; see docs_site/ for that
├── docs_site/                # This documentation site's Markdown source (mkdocs.yml → docs_dir)
├── tools/                    # Repository tools: golden-run gate for solver changes (tools/golden/)
├── .github/workflows/ci.yml  # GitHub Actions: fast suite on every push/PR, medium suite on main
├── runs/                     # Per-run I/O directories (gitignored)
├── build/                    # Compiled binaries (gitignored)
├── logs/                     # SLURM stdout/stderr logs (gitignored)
├── Makefile                  # make build / build-video / run / submit / clean
├── mkdocs.yml                # This documentation site's configuration
└── pyproject.toml            # Python dependencies (managed by uv)
```
