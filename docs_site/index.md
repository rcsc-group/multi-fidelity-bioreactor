# Rocking Bioreactor 2D — Simulation Suite

Two-phase CFD solver for a rocking bioreactor, implemented in [Basilisk](http://basilisk.fr/).
Developed at the Harris Lab (Brown University) in collaboration with the Cimpeanu group (Warwick).

Publication: [doi: 10.1016/j.ijmultiphaseflow.2025.105375](https://www.sciencedirect.com/science/article/pii/S0301932225002538) | preprint: [arXiv: 2504.05421](https://arxiv.org/abs/2504.05421)

The solver resolves two-phase (VOF) hydrodynamics and dissolved-oxygen transport
(Henry's law), producing kLa (volumetric oxygen mass-transfer coefficient) and
shear-stress KPIs as optimization objectives. A multi-fidelity Bayesian
optimization suite sits on top, using a KRR-LR-GPR surrogate to trade off
cheap low-fidelity screening runs against expensive high-fidelity corrections.

## Where to start

- New to the project? Read the [Glossary](glossary.md) first.
- Setting up an environment? See [Setup](setup.md).
- Want to run something? Pick a [workflow](workflows/single-run.md) that matches
  what you're trying to do — a single run, a parameter sweep, or the full BO loop.
- Looking up a field or file format? See the [Reference](reference/params.md) section.

## Repository layout at a glance

```
├── src/         # Basilisk solver (BioReactor.c + headers)
├── scripts/     # Python orchestration (simulate, sweep, BO loop, postprocess, ...)
├── config/      # YAML/JSON configs for each workflow
├── tests/       # Unit, integration, and numerical-verification tests
├── experiments/ # Sweep metadata, chain manifests, generated figures
├── docs/        # Reference papers (Kim et al. 2024/2025) and the canonical validation case
└── mkdocs.yml   # This documentation site (see Reference → Project structure for the full tree)
```

See [Project structure](reference/project-structure.md) for the complete annotated tree.
