# Rocking Bioreactor 2D — Simulation Suite

[![CI](https://github.com/rcsc-group/multi-fidelity-bioreactor/actions/workflows/ci.yml/badge.svg)](https://github.com/rcsc-group/multi-fidelity-bioreactor/actions/workflows/ci.yml)

Two-phase CFD solver for a rocking bioreactor, implemented in [Basilisk](http://basilisk.fr/).
Developed at the Harris Lab (Brown University) in collaboration with the Cimpeanu group (Warwick).
Publication: [doi: 10.1016/j.ijmultiphaseflow.2025.105375](https://www.sciencedirect.com/science/article/pii/S0301932225002538) | preprint: [arXiv: 2504.05421](https://arxiv.org/abs/2504.05421)

The solver resolves two-phase (VOF) hydrodynamics and dissolved-oxygen transport
(Henry's law), producing kLa (volumetric oxygen mass-transfer coefficient) and
shear-stress KPIs as optimization objectives. A multi-fidelity Bayesian
optimization suite sits on top — a KRR-LR-GPR surrogate trained on cheap
low-fidelity screening runs plus a smaller number of expensive high-fidelity
corrections, with an Expected Improvement acquisition function choosing where
to sample next.

## Quickstart

```bash
pip install uv && uv sync   # Python dependencies
make build                  # compiles build/BioReactor (needs Basilisk's qcc on PATH)
uv run python -m pytest tests/   # fast unit suite, ~20s, no binary or SLURM needed
```

## Documentation

Full setup, workflow tutorials, and reference material live in [`docs_site/`](docs_site/index.md)
(built with [MkDocs](https://www.mkdocs.org/); run `uv run --group docs mkdocs serve` for a live local preview):

- [Glossary](docs_site/glossary.md) — start here if you're new to CFD/HPC terminology used throughout
- [Setup](docs_site/setup.md)
- Workflows: [single run](docs_site/workflows/single-run.md) ·
  [chained sweep](docs_site/workflows/chained-sweep.md) ·
  [batch sampling](docs_site/workflows/batch-sampling.md) ·
  [Bayesian optimization](docs_site/workflows/bayesian-optimization.md) ·
  [JSON multi-param sweep](docs_site/workflows/json-sweep.md)
- Reference: [params.json](docs_site/reference/params.md) ·
  [output files](docs_site/reference/output-files.md) ·
  [fidelity guide](docs_site/reference/fidelity-guide.md) ·
  [project structure](docs_site/reference/project-structure.md)
- [Test suite](docs_site/testing.md)

## References

- Kim M., Harris D.M., Cimpeanu R. (2025). *Modelling of oxygen transfer in rocking bioreactors.* Int. J. Multiphase Flow. [doi: 10.1016/j.ijmultiphaseflow.2025.105375](https://www.sciencedirect.com/science/article/pii/S0301932225002538)
- Basilisk CFD: http://basilisk.fr
