# Rocking Bioreactor 2D — Simulation Suite

Two-phase CFD solver for a rocking bioreactor, implemented in [Basilisk](http://basilisk.fr/).
Developed at the Harris Lab (Brown University) in collaboration with the Cimpeanu group (Warwick).

Publication: [doi: 10.1016/j.ijmultiphaseflow.2025.105375](https://www.sciencedirect.com/science/article/pii/S0301932225002538) | preprint: [arXiv: 2504.05421](https://arxiv.org/abs/2504.05421)

The solver resolves two-phase (VOF) hydrodynamics and dissolved-oxygen transport
(Henry's law), producing kLa (volumetric oxygen mass-transfer coefficient) and
shear-stress KPIs as optimization objectives. A multi-fidelity Bayesian
optimization suite sits on top, using a KRR-LR-GPR surrogate to trade off
cheap low-fidelity screening runs against expensive high-fidelity corrections.

## How this documentation is organized

- **[Tutorials](tutorials/first-simulation.md)** — learning by doing, start to
  finish, with real output shown. Start here if you're new to the project.
- **How-to guides** — task-focused instructions assuming you already know the
  basics: [sweep one parameter](how-to/sweep-one-parameter.md),
  [sweep any parameter combination](how-to/sweep-json-multi-param.md),
  [batch sampling](how-to/batch-sampling.md),
  [run the BO loop](how-to/run-bo-loop.md),
  [generate videos](how-to/generate-videos.md),
  [choose a fidelity level](how-to/choose-fidelity.md),
  [diagnose a stalled chain](how-to/diagnose-stalled-chain.md).
- **[Reference](reference/params.md)** — dry, structured lookup material:
  params.json fields, output files, the fidelity table, scripts, project layout.
- **Explanation** — the "why": [checkpoint restart and warm-start chains](explanation/checkpoint-restart.md),
  [non-dimensionalization](explanation/non-dimensionalization.md),
  [multi-fidelity Bayesian optimization](explanation/multi-fidelity-bo.md),
  and an honest account of [where validation against Kim et al. (2024) actually stands](explanation/kim-et-al-validation.md).

New to CFD/HPC terminology used throughout? Read the [Glossary](glossary.md) first.
Setting up an environment? See [Setup](setup.md).

## Repository layout at a glance

```
├── src/         # Basilisk solver (BioReactor.c + headers)
├── scripts/     # Python orchestration (simulate, sweep, BO loop, postprocess, ...)
├── examples/    # Runnable, self-contained tutorial companions
├── config/      # YAML/JSON configs for each workflow
├── tests/       # Unit, integration, and numerical-verification tests
├── experiments/ # Sweep metadata, chain manifests, generated figures
├── docs/        # Reference papers (Kim et al. 2024/2025) and the canonical validation case
└── mkdocs.yml   # This documentation site (see Reference → Project structure for the full tree)
```

See [Project structure](reference/project-structure.md) for the complete annotated tree.
