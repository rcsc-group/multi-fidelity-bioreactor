# Rocking Bioreactor 2D 

[![CI](https://github.com/rcsc-group/multi-fidelity-bioreactor/actions/workflows/ci.yml/badge.svg)](https://github.com/rcsc-group/multi-fidelity-bioreactor/actions/workflows/ci.yml)
[![Docs](https://github.com/rcsc-group/multi-fidelity-bioreactor/actions/workflows/docs.yml/badge.svg)](https://rcsc-group.github.io/multi-fidelity-bioreactor/)

**📖 Documentation: [rcsc-group.github.io/multi-fidelity-bioreactor](https://rcsc-group.github.io/multi-fidelity-bioreactor/)**

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

Full setup, workflow tutorials, and reference material live on the
[documentation site](https://rcsc-group.github.io/multi-fidelity-bioreactor/)
(source in [`docs_site/`](docs_site/index.md), built with [MkDocs](https://www.mkdocs.org/);
run `uv run --group docs mkdocs serve` for a live local preview — auto-deployed
to Pages on every push to `main`):

- [Glossary](https://rcsc-group.github.io/multi-fidelity-bioreactor/glossary/) — start here if you're new to CFD/HPC terminology used throughout
- [Setup](https://rcsc-group.github.io/multi-fidelity-bioreactor/setup/)
- Tutorials: [your first simulation](https://rcsc-group.github.io/multi-fidelity-bioreactor/tutorials/first-simulation/) ·
  [your first sweep](https://rcsc-group.github.io/multi-fidelity-bioreactor/tutorials/first-sweep/) ·
  [your first optimization loop](https://rcsc-group.github.io/multi-fidelity-bioreactor/tutorials/first-optimization-loop/) (no Basilisk build needed)
- How-to guides: [sweep one parameter](https://rcsc-group.github.io/multi-fidelity-bioreactor/how-to/sweep-one-parameter/) ·
  [sweep any parameter combination](https://rcsc-group.github.io/multi-fidelity-bioreactor/how-to/sweep-json-multi-param/) ·
  [run a batch DoE](https://rcsc-group.github.io/multi-fidelity-bioreactor/how-to/batch-sampling/) ·
  [run/resume the BO loop](https://rcsc-group.github.io/multi-fidelity-bioreactor/how-to/run-bo-loop/) ·
  [diagnose a stalled chain](https://rcsc-group.github.io/multi-fidelity-bioreactor/how-to/diagnose-stalled-chain/)
- Explanation: [checkpoint restart](https://rcsc-group.github.io/multi-fidelity-bioreactor/explanation/checkpoint-restart/) ·
  [non-dimensionalization](https://rcsc-group.github.io/multi-fidelity-bioreactor/explanation/non-dimensionalization/) ·
  [multi-fidelity BO](https://rcsc-group.github.io/multi-fidelity-bioreactor/explanation/multi-fidelity-bo/) ·
  [validating against Kim et al. (2024)](https://rcsc-group.github.io/multi-fidelity-bioreactor/explanation/kim-et-al-validation/)
- Reference: [params.json](https://rcsc-group.github.io/multi-fidelity-bioreactor/reference/params/) ·
  [output files](https://rcsc-group.github.io/multi-fidelity-bioreactor/reference/output-files/) ·
  [fidelity guide](https://rcsc-group.github.io/multi-fidelity-bioreactor/reference/fidelity-guide/) ·
  [scripts](https://rcsc-group.github.io/multi-fidelity-bioreactor/reference/scripts/) ·
  [project structure](https://rcsc-group.github.io/multi-fidelity-bioreactor/reference/project-structure/)
- [Test suite](https://rcsc-group.github.io/multi-fidelity-bioreactor/testing/)

## References

- Kim M., Harris D.M., Cimpeanu R. (2025). *Modelling of oxygen transfer in rocking bioreactors.* Int. J. Multiphase Flow. [doi: 10.1016/j.ijmultiphaseflow.2025.105375](https://www.sciencedirect.com/science/article/pii/S0301932225002538)
- Basilisk CFD: http://basilisk.fr
