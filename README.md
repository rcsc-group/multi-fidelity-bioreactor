# Rocking bioreactor: two-phase CFD and multi-fidelity optimisation

[![CI](https://github.com/rcsc-group/multi-fidelity-bioreactor/actions/workflows/ci.yml/badge.svg)](https://github.com/rcsc-group/multi-fidelity-bioreactor/actions/workflows/ci.yml)
[![Docs](https://github.com/rcsc-group/multi-fidelity-bioreactor/actions/workflows/docs.yml/badge.svg)](https://rcsc-group.github.io/multi-fidelity-bioreactor/)

<video src="docs_site/assets/img/hero-rocking-l10-lab.mp4" controls muted loop playsinline width="600">
  Rocking bioreactor simulation at level 10 (1024x1024 cells), 7 degrees, 32.5 rpm, lab frame, coloured by speed.
</video>

A rocking bag of cell-culture medium: water and air, rocked about an axis. This repository has

- a **solver** ([Basilisk](http://basilisk.fr/), two-phase VOF with an embedded wall) that computes the
  flow, the oxygen transfer (Henry's law) and a passive tracer. It reports four quantities: shear stress,
  oxygen transfer rate kLa, energy dissipation rate and mixing time;
- a **multi-fidelity surrogate and Bayesian optimisation** layer: many cheap coarse-grid runs plus a few
  expensive fine-grid runs, combined into one prediction of those quantities over the rocking angle and
  rate.

The solver started as the driver code of Kim, Harris and Cimpeanu
([rcsc-group/BioReactor](https://github.com/rcsc-group/BioReactor)) and has since been extended (any grid
level, checkpoint restart, more outputs, flux record/replay).

**Documentation:** <https://rcsc-group.github.io/multi-fidelity-bioreactor/> (source in `docs_site/`).

## Quick start

```bash
pip install uv && uv sync              # Python dependencies
make build                             # build/BioReactor (needs Basilisk's qcc)
uv run python -m pytest tests/         # fast unit tests, no binary or SLURM needed
```

Then follow [Your first simulation](https://rcsc-group.github.io/multi-fidelity-bioreactor/tutorials/first-simulation/).
On a shared cluster, run simulations through SLURM (`make submit PARAMS=...`), never on a login node.

## Where do I start?

| I want to... | Start at |
|---|---|
| understand the physics and the method | [Explanation pages](https://rcsc-group.github.io/multi-fidelity-bioreactor/explanation/non-dimensionalization/) |
| run one simulation | [`docs_site/tutorials/first-simulation.md`](docs_site/tutorials/first-simulation.md), `config/` |
| change the solver | [`src/README.md`](src/README.md), then [solver pitfalls](docs_site/developer/solver-pitfalls.md) |
| run a sweep or a campaign | [`docs_site/how-to/`](docs_site/how-to/), `scripts/README.md` |
| fit the surrogate / run the optimiser | [`docs_site/explanation/multi-fidelity-bo.md`](docs_site/explanation/multi-fidelity-bo.md), `scripts/bo2d/` |
| rebuild a figure | `experiments/README.md` |
| run the tests | [`docs_site/testing.md`](docs_site/testing.md) |

## Repository map

| Folder | Contents |
|---|---|
| `src/` | The Basilisk solver: `BioReactor.c` plus one header per concern. Map and event order: `src/README.md`. |
| `scripts/` | Python: run and chain simulations, postprocess, campaigns (`campaign2d/`), surrogate and optimisation (`bo2d/`, `mfbml_local/`), figures. |
| `config/` | Parameter spaces, sweep and chain configurations, SLURM templates. |
| `tests/` | Unit tests (default), `medium`/`hpc` tests (CI or cluster only), `verification/` physics checks. |
| `tools/` | Repository tools: golden-run gate for solver changes (`tools/golden/`), repository audit. |
| `experiments/` | Results that figures and the diary use (JSON/CSV), and figures. |
| `docs_site/` | Public documentation (MkDocs). |
| `docs/` | Internal notes: specs, literature, meetings, solver history. Not published. |
| `talks/` | Talk material. |
| `examples/`, `sandbox/` | A synthetic optimisation demo; a small MPI restart reproducer. |
| `diary.md` | Lab notebook: what each experiment found (older entries: `diary_archive_to_2026-10-07.md`). |
| `BACKLOG.md` | Open work. |
| `runs` | Link to the run data on the OSCAR cluster (not portable). |

## Reference

- Kim M., Harris D.M., Cimpeanu R. (2025). *Modelling of oxygen transfer in rocking bioreactors.*
  Int. J. Multiphase Flow. [doi:10.1016/j.ijmultiphaseflow.2025.105375](https://www.sciencedirect.com/science/article/pii/S0301932225002538)
- Basilisk: <http://basilisk.fr>
