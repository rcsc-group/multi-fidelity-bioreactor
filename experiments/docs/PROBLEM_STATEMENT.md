# Rocking Bag Bioreactor — Operating Condition Design Problem

## Problem

Single-use rocking bag bioreactors are limited by oxygen transfer. The platform
rocks sinusoidally, driving a slow time-averaged secondary circulation — **steady
streaming** — that transports dissolved oxygen from the headspace into the liquid
and mixes the culture. The oxygen transfer rate is characterised by the
**volumetric oxygen transfer coefficient** *k*L*a* [h⁻¹]:

```
dC*(t)/dt = kLa × (1 − C*(t))
```

where *C*\*(t) ∈ [0, 1] is dimensionless dissolved oxygen saturation.

The design question is: **what platform motion and bag geometry maximise *k*L*a*
at the 25% saturation crossing?**

Kim et al. (2025) — the current state of the art — explored only two parameters
(rocking frequency and amplitude) at a single geometry and fill level. They
report a best of **267 h⁻¹** at 37.5 RPM, θ = 7°. Fourteen dimensions of the
operating space were left untouched: fill level, bag shape, multi-harmonic
waveforms, and horizontal oscillation. This benchmark targets those dimensions.

---

## Evaluator

The evaluator is the compiled `BioReactor` Basilisk direct Navier–Stokes solver
(2D VOF, embedded boundary, quadtree AMR) wrapped by `scripts/postprocess.py`.
It is exposed through the f3dasm interface via `BioreactorDataGenerator` in
`workspace/`:

```python
from f3dasm import ExperimentData
from data_generator import BioreactorDataGenerator

data = ExperimentData(domain=domain)
data.sample(sampler="latin_sampler", n_samples=N)
data.run(data_generator=BioreactorDataGenerator())
```

Each evaluation writes a `params.json`, runs the solver, and extracts twelve
KPIs from `tr_oxy.dat`, `vol_frac_interf.dat`, and `normf.dat`. Simulation cost
scales as N⁴ where N = 2^fidelity. A fidelity-7 run takes ~4 h on 16 MPI ranks
on OSCAR. **`fidelity` is a solver parameter, not a design variable.**

The full evaluation stack — DoE, surrogate training, acquisition, SLURM
submission, checkpoint restart, and health monitoring — is implemented in:

| Script | Purpose |
|---|---|
| `scripts/sample.py` | Space-filling DoE batch (Latin hypercube, Sobol, random) |
| `scripts/sweep.py` | Multi-parameter JSON sweep with automatic checkpoint grouping |
| `scripts/chain.py` | Checkpoint-restart chain for sequential designs sharing geometry |
| `scripts/loop.py` | Multi-fidelity Bayesian optimisation loop (KRR-LR-GPR surrogate) |
| `scripts/train_surrogate.py` | Fit KRR-LR-GPR from ExperimentData store |
| `scripts/suggest.py` | Acquisition function (EI, UCB) evaluated on the surrogate |
| `scripts/postprocess.py` | Extract all 12 KPIs from a completed run |
| `scripts/health_report.py` | Sweep-level health check (NaN counts, dtmix ordering) |

---

## Design space

### Notation

The platform motion is:

```
θ(t) = Σ_{k=0}^{n_harmonics-1}  theta_max[k] · sin((k+1)·omega_b·t + phi_angular[k])
x(t) = Σ_{k=0}^{n_harmonics-1}  amplitude_h[k] · sin((k+1)·omega_h·t + phi_horizontal[k])
```

`phi_angular[0] = 0` always (time-origin reference; not a free parameter).
The bag cross-section is `|x/a|^n + |y/b|^n ≤ 1`.

Inactive harmonic components are padded to zero and validated by
`postprocess.validate_params()` against `config/param_space.yaml`.

### Free parameters (18 total at n_harmonics = 3)

| Parameter | Bounds | Unit | Notes |
|---|---|---|---|
| `omega_b` | [1.57, 6.28] | rad/s | 15–60 RPM |
| `n_harmonics` | {1, 2, 3} | — | integer; gates active components |
| `theta_max[0..2]` | [2.0, 7.0] each | deg | active harmonics only; inactive = 0 |
| `phi_angular[1..2]` | [0, 2π] each | rad | waveform shape; only active when n_harmonics ≥ 2/3 |
| `omega_h` | [0.0, 6.28] | rad/s | 0 = pure rocking |
| `amplitude_h[0..2]` | [0.0, 0.05] each | m | active when omega_h > 0 |
| `phi_horizontal[0..2]` | [0, 2π] each | rad | phase of horizontal vs. rocking |
| `geometry.a` | [0.15, 0.35] | m | bag semi-length |
| `geometry.b` | [0.05, 0.15] | m | bag semi-height; must satisfy b < a |
| `geometry.n` | [2.0, 8.0] | — | 2 = ellipse, ≥ 6 ≈ rectangle |
| `fill_level` | [0.3, 0.7] | — | liquid volume fraction |

**Effective dimensionality by subspace:**

| Subspace | Active parameters | Dim |
|---|---|---|
| Kim et al. (frequency + amplitude only) | `omega_b`, `theta_max[0]` | 2 |
| Single-harmonic, default geometry | + `fill_level`, `geometry.{a,b,n}` | 6 |
| Multi-harmonic, no horizontal | + `theta_max[1..2]`, `phi_angular[1..2]`, `n_harmonics` | 11 |
| Full space | all above + `omega_h`, `amplitude_h[0..2]`, `phi_horizontal[0..2]` | 18 |

---

## Outputs

`postprocess.py` returns twelve quantities per completed run:

| Key | Description | Unit |
|---|---|---|
| `kLa_10` | O₂ transfer rate at C\*=10%, 5-pt log-linear fit | 1/t_nd |
| `kLa_25` | O₂ transfer rate at C\*=25% ← **primary objective** | 1/t_nd |
| `kLa_50` | O₂ transfer rate at C\*=50% | 1/t_nd |
| `kLa_inst_10` | Instantaneous kLa at C\*=10% | 1/t_nd |
| `kLa_inst_25` | Instantaneous kLa at C\*=25% | 1/t_nd |
| `kLa_inst_50` | Instantaneous kLa at C\*=50% | 1/t_nd |
| `dtmix_0.50` | Time to 50% passive-tracer mixing | seconds |
| `dtmix_0.75` | Time to 75% mixing | seconds |
| `dtmix_0.95` | Time to 95% mixing | seconds |
| `vor_mean` | Period-averaged mean absolute vorticity (streaming proxy) | 1/s |
| `vel_rms_qss` | RMS velocity in quasi-steady window | non-dim |
| `kla_fit_rmse_25` | RMSE of kLa_25 log-linear fit | — |

`vel_rms_qss` is the convergence diagnostic from Kim et al. Appendix A.
`kla_fit_rmse_25 < 0.005` indicates a reliable kLa_25 estimate.
`vor_mean` is a cheap low-fidelity proxy for kLa that can guide early exploration
before oxygen curves are meaningful.

### Unit conversion

The solver runs in non-dimensional time. The dimensional kLa is recovered as:

```
T_per  = 2π / omega_b
V_char = (geometry.a / 4) × (geometry.b  +  0.5 × geometry.a × tan(theta_max[0]))
U_bio  = V_char / (geometry.b / 2) / T_per
T_bio  = geometry.a / U_bio          ← characteristic time; NOT equal to T_per

kLa_25 [1/s]  =  kLa_25_nd / T_bio
kLa_25 [h⁻¹] =  kLa_25 [1/s] × 3600
```

This conversion is implemented in `postprocess._t_scales()`. Note:
`T_bio ≠ T_per`. For the default geometry (a=0.25, b=0.071, θ=7°),
T_bio/T_per = 1.645 — a constant geometric ratio. The incorrect formula
`kLa_nd × omega_b / (2π)` overstates kLa_dim by this factor and must not
be used. `BioreactorDataGenerator` returns `kLa_25_dim` in h⁻¹ using
the correct T_bio conversion.

---

## Objective

```
maximise   kLa_25_dim [h⁻¹]   subject to  no overflow
```

subject to `kla_fit_rmse_25 < 0.005` (discard noisy estimates).

---

## Constraint

**No overflow.** `posY_max < geometry.b / geometry.a` must hold at all times.
`BioreactorDataGenerator` returns `kLa_25_dim = NaN` and `overflow = True`
for constraint-violating designs.

---

## Reference baselines (fidelity 7, this evaluator)

All comparisons are made within this evaluator at fixed fidelity to be
internally consistent.

| Label | Condition | kLa_25_dim |
|---|---|---|
| **B-Kim** | Kim et al. operating point: `omega_b=3.927`, `theta_max=[7,0,0]`, `fill=0.5`, `geometry={a:0.25, b:0.071, n:8}` | **38.2 h⁻¹** |
| **B-2D** | Best from axis-aligned 2D sweeps (fill × ω and theta × ω, both at `geometry={a:0.25, b:0.071, n:8}`): `omega_b=1.833`, `theta_max=[7,0,0]`, `fill=0.3` | **88.6 h⁻¹** |

B-Kim is the anchor: it establishes where the standard operating condition sits
within this evaluator and is comparable in methodology (same fidelity, same
solver) to any new design proposed here.

B-2D is the current ceiling reachable by searching only two parameters at a
time. Any strategy that cannot beat B-2D is not using the high-dimensional
space.

---

## Multi-fidelity strategy

The recommended search strategy uses all four fidelity tiers:

```
Phase 1 — Landscape (fidelity 4–5, ~seconds–minutes each)
    scripts/sample.py  config/sample_config.yaml
    Large LHS DoE over all 18 dimensions.
    Use vor_mean as proxy objective; kLa_nd too noisy at this resolution.
    Target: ~500–1000 evaluations mapping the full space.

Phase 2 — Surrogate + acquisition (fidelity 6, ~10 min each)
    scripts/train_surrogate.py  →  scripts/suggest.py  →  scripts/simulate.py
    Fit KRR-LR-GPR multi-fidelity surrogate on fidelity 4+5 data.
    Acquire candidate designs via EI/UCB; evaluate at fidelity 6.
    Target: ~50–100 evaluations refining promising regions.

Phase 3 — Confirmation (fidelity 7, ~4 h each, 16 MPI)
    scripts/sweep.py  config/sweep_<candidate>.json
    Confirm the top-5 surrogate candidates at fidelity 7.
    Chain restart (scripts/chain.py) where multiple candidates share geometry.
    Target: ~10–20 evaluations; these are the reportable results.

Phase 4 — Mechanism (fidelity 7, targeted)
    Axis-aligned sweeps around the winning design to isolate which parameters
    drive the gain (sweep one axis, fix all others at optimum).
    Required for physical interpretation.
```

`loop.py` automates Phases 2–3 end-to-end when run with a suitable
`bo_config.yaml`. Health checks after each phase via `scripts/health_report.py`.

---

## Success criterion

A solution to this problem must deliver:

1. **A specific design vector** — all 18 parameters fully specified, valid under
   `param_space.yaml`, confirmed at fidelity 7 by this evaluator.

2. **kLa_25_dim > 88.6 h⁻¹** (B-2D) with `kla_fit_rmse_25 < 0.005` and
   no overflow.

3. **Improvement from at least one axis Kim et al. did not explore.** The
   winning design must differ from B-Kim and B-2D in at least one of:
   `geometry.n ≠ 8`, `fill_level ≠ 0.5`, `n_harmonics ≥ 2`, `omega_h > 0`,
   or `geometry.{a,b}` off-default. A design that only refines frequency and
   amplitude is not a valid solution.

4. **Mechanism attribution.** For the winning design, the following must be
   reported from KPI data:
   - `vor_mean` relative to B-Kim and B-2D (quantifies streaming enhancement)
   - `dtmix_0.95` (confirms culture homogeneity)
   - Which single parameter change from B-2D produces the largest marginal gain
     in kLa_25_dim (identified by a targeted axis sweep at fidelity 7)

5. **Evaluation budget accounting.** Total evaluations at each fidelity tier
   and the search strategy used.

---

## Physical context

The five axes left unexplored by Kim et al. and their expected mechanisms:

**Fill level.** Lower fill → higher headspace-to-liquid ratio → more oxygen
interface per unit liquid volume. Also confines the streaming vortex pair into a
shallower layer, increasing local shear. Preliminary sweep (`sweep_fb_fill_l7_v2`)
shows fill = 0.3 outperforms fill = 0.5 at all tested frequencies (88.6 vs
38.2 h⁻¹ at matched omega_b).

**Multi-harmonic waveform.** A second rocking harmonic at 2ω_b with phase φ₁
breaks the time-symmetry of the forcing. Via quadratic nonlinearity in the
Navier–Stokes equations, asymmetric oscillatory forcing generates a stronger
time-averaged streaming component — directly increasing vor_mean and kLa.
Whether this effect survives at finite Reynolds number is untested.

**Horizontal oscillation.** Combining angular rocking (θ) with horizontal
translation (x) produces a Lissajous platform trajectory. Specific phase
relationships between the two motions break the left–right vortex symmetry
and can drive net directed recirculation — a qualitatively different transport
mode from the symmetric counter-rotating vortex pair produced by pure rocking.

**Bag cross-section shape.** Near-rectangular bags (geometry.n ≥ 6) generate
sharper corner vortices than elliptical ones (geometry.n = 2). Whether corner
vorticity enhances or disrupts the streaming pair is unknown; it is likely
geometry-dependent and interacts with fill level.

**Bag aspect ratio.** Varying geometry.a and geometry.b at fixed fill level
changes both the absolute scale and the aspect ratio of the liquid layer.
Thinner, longer bags (high a/b) may favour wave propagation over the surface;
shorter, deeper bags may favour bulk recirculation. Scale also sets the
Womersley number Wo = (geometry.b) × sqrt(omega_b / nu), which controls the
oscillatory boundary layer thickness relative to the liquid depth.
