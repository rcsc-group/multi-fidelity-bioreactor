# Multi-fidelity prediction of a grid-converged quantity: problem definition and algorithm

Draft 7, 2026-10-03. It replaces draft 1 (chat, 2026-10-02) and drafts 2–6. The adversarial reviewer rejected drafts 2–5 in rounds 1–4 (1 fatal and 6 major holes, then 3, 2 and 2 major holes) and approved draft 6 in round 5 with 6 minor holes. Draft 7 closes those minor holes.

**Status labels:**
- **[lit]**: I read it in the cited paper in this session.
- **[data]**: our measurement; the script is named.
- **[inference]**: my reasoning; it is not tested.
- **[assumption]**: a choice; it is not tested.

---

## 0. Main findings that shape this draft

| # | Finding | Evidence |
|---|---|---|
| F1 | **Provisional: the grid-converged mixing time is not identifiable from L6–L9.** For three consecutive levels the increment ratio \(R = (z_2 - z_1)/(z_3 - z_2)\) must be about \(2^p > 1\) for convergence (see the note below the table). This uses per-rpm triplets. The pooled G0 test of Section 4 is first action 1 and has not run yet. | [data] Fig. D; `scripts/diag_observed_order.py` prints R and \(r_\infty\) |
| F2 | For kLa (no transform), convergence is **plausible, not shown**. L6–L8 is monotone at 90–100% of rpm (C* = 25, 50%). L7–L9 is monotone at only 56–60%, and 4 of 10 rpm are oscillatory. So convergence depends on the QoI and on the transform, and the method must test it for each QoI. | [data] `observed_order_kla.png` |
| F3 | "Kink at χ ≈ 0.75" (draft 1) is **falsified**. The late/early decay-rate ratio is 0.36–3.65, and its direction changes with rpm and with level. | [data] Fig. A |
| F4 | Stopping at χ = 0.5 saves only 1.2–1.5× (L8) and 1.5–2.2× (L9), not 3–9×. The cause is the 80-cycle spin-up, which is 40–55% of an L9 run. | [data] diary 2026-10-03 |
| F5 | The cost per level is steep. Per simulated second: ×7, ×9, ×21. Per probe to \(\Delta t_{0.95}\): ×17, ×24, ×50. The second set is larger because \(\Delta t\) also grows. | [data] Fig. C |

Details of F1:
- \(z = \Delta t\): R = 0.2–0.4 (median), so the increments grow.
- \(z = \log \Delta t\): R = 0.6–1.2, so the increments are constant. This is the transform that predicts the next level best (Fig. B).
- \(z = 1/\Delta t\): R = 1.5–3.4, monotone at 70–89% of rpm. But the Richardson limit \(r_\infty / r_{L9}\) from L7–L9 is −2.3 to 0.98 across rpm (medians 0.05, 0.30 and 0.46 at χ = 0.50, 0.75 and 0.95). 9 of 22 values are negative, which is unphysical. So the limit is compatible with \(r_\infty = 0\), which means an infinite mixing time.
- The pooled rate fit with one shared \(p\) (diary 2026-10-02) identifies \(p = 1.41 \pm 0.31\) (χ = 0.95), but its \(r_\infty\) is small against its uncertainty and negative at some rpm. So pooling identifies \(p\) but not an admissible limit.
- **[inference]** One possible cause is physical. With \(D = 0.44 \times 10^{-9}\) m²/s (`src/BioReactor.c:201`) and \(L = 0.25\) m, the Péclet number is about \(10^7\). The Batchelor scale is then about 50 µm, against an L9 cell of 0.49 mm. On every level we run, numerical diffusion may set the final mixing. If so, the asymptotic range starts near L12–L13, which we cannot afford (F5).

![Fig. D](../experiments/multifidelity/observed_order_dtmix.png)

*Fig. D. Increment ratio R per rpm for L6–L8 (circles) and L7–L9 (squares). Rows: \(\Delta t\), \(\log \Delta t\), \(1/\Delta t\). Dashed: R = 1, no convergence. Dotted: R = 2.8, order 1.5.*

![Fig. A](../experiments/multifidelity/chi_decay_rate_shape.png)

*Fig. A. The period-averaged decay rate \(\lambda(\chi)/\lambda(0.5)\), with \(\lambda = -d\ln(1-\chi)/dt\), for L6–L9.*

**Consequence:** the method must first decide whether its target is identifiable (Gate G0, Section 4). Section 8 shows that, for mixing time at \(h = 0\) with L6–L9:
- two of the three transforms fail G0;
- the third (the rate) can pass G0, but its limit is not bounded away from zero, so the Step 4 forecast is expected to report infeasibility.

Both statements are provisional until the pooled fit (action 1).

---

## 1. Setting and notation

- \(x \in X \subset R^d\): the QoI coordinates. \(\Sigma \subset X\) is the region of interest. \(\Sigma_N\) is a scrambled Sobol set of \(100 \cdot d\) points in \(\Sigma\). All "max over \(\Sigma\)" values are computed on \(\Sigma_N\) [assumption].
- \(h\): the cell size. The levels are \(h_\ell = h_0 2^{-\ell}\). \(h_{min}\) is the finest level that has data. \(h_c\) is the coarsest level of any candidate level set (fixed for all structures).
- **One normalisation:** inside every kernel and prior, \(h\) enters as \(\bar h = h/h_c \in (0, 1]\). Figures use \(h/h_{ref}\) with \(h_{ref}\) the L10 cell, for display only.
- All coordinates of \(x\) are scaled to [0, 1] over \(\Sigma\) inside the kernels.
- **Target resolution \(h^\star\).** The default is \(h^\star = 0\), the grid-converged value. The method also supports \(h^\star > 0\), for example the resolution of a reference study. Then the target is \(f(x, h^\star)\).
- A **probe** is one run, \(a = (u, h, T)\): controls \(u\), level \(h\), run length \(T\).
  - It returns a vector \(y_a\) at the coordinate set \(O(a) \subset X\).
  - Generic case: \(O(a) = \{u\}\).
  - Nested case: \(x_k = (\chi_k, u)\) for each grid value \(\chi_k\) reached before \(T\). So \(O(a)\) depends on the outputs.
- \(g\): the output transform. \(c(a)\): the cost in core-hours. \(C\): the budget.
- \(\varepsilon\): the tolerance, on the **physical** scale. It is either relative (\(\varepsilon_{rel}\)) or absolute (\(\varepsilon_{abs}\)).

**Scope.** The target \(f(x, 0)\) is the converged output of *the simulation model*. Model-form error, such as 2D instead of 3D or the turbulence treatment, is outside this method. To include it, the method needs experimental data and a model-discrepancy term (question Q5).

---

## 2. Model. The structure is generic; only the priors in Section 8 are problem-specific.

### 2.1 Observation model

For output \(k\) of probe \(a\):

$$ g(y_{a,k}) = \mu(x_k) + \delta(x_k, h_a) + e_{a,k} $$

- \(\mu\) is the converged mean in \(g\)-space.
- \(\delta\) is the discretisation error, with \(\delta(x, 0) = 0\).
- \(e_a\) is the noise vector of the run.
- **[lit]** This is the Tuo–Wu–Yu decomposition, restated in CONFIG 2209.13748 eq 6–8; Stroh 1605.02561 eq 1–2 adds noise.
- All levels enter one likelihood. No level is the truth.

### 2.2 Converged mean

$$ \mu \sim GP(\phi(x)^T \beta, \; \sigma_\mu^2 k_\mu(x, x')) $$

- \(\beta\) has a flat prior (universal kriging) **only when there is no link**. Then the likelihood is Gaussian in \(\beta\), and the posterior is proper when the basis matrix \(F\) (rows \(\phi(x_i)^T\) at the data sites) has full column rank.
  - A basis term that does not vary in the data, for example \(\log \theta_{max}\) when all runs have one angle, is collinear with the intercept. Its coefficient then gets a proper \(N(0, 1)\) prior until the term varies at 2 or more levels.
- \(k_\mu\) is a Matérn-5/2 kernel with one length scale per input dimension (ARD) [assumption]. The same choice is used for \(k_x\) below.
- **Known admissible range goes into the prior, not into a gate.** If the converged value must lie in a known range, write \(\mu = \psi(\tilde\mu)\) with a link \(\psi\) and a GP on \(\tilde\mu\). For example, for \(g = -1/\Delta t\) the converged rate must be positive, so \(\mu = -\exp(\tilde\mu)\). \(\delta\) stays additive in \(g\)-space. The likelihood is then not Gaussian in \(\tilde\mu\), and Section 2.5 samples it.
  - **With a link, \(\beta\) must have a proper prior.** With a flat \(\beta\), the intercept can go to \(-\infty\): then \(\mu \to 0\) everywhere, the likelihood tends to a positive constant, and the posterior is improper. We use \(\beta \sim N(b_{phys}, \mathrm{diag}(s_b^2))\), with \(b_{phys}\) from a physical estimate (Section 8).
  - The link constrains only \(\mu\). For a finite target \(h^\star > 0\), \(f(x, h^\star) = \psi(\tilde\mu) + \delta(x, h^\star)\) can still leave the range. The method reports the probability that \(f(x, h^\star)\) is outside the range, and flags values above 5%.

### 2.3 Discretisation error

$$ \delta \sim GP(0, \; \sigma_\delta^2 k_x(x, x') k_h(h, h')) $$

There are two candidate families for \(k_h\). Section 2.6 weights them.

- **TWY2 (Richardson type):** \(k_h = (h h')^p c_\nu(h - h'; \ell_h)\), with \(c_\nu\) a Matérn correlation.
  - [lit] Bect et al. 2103.14559 §3 and Prop. 3: \(\delta(h) = A h^p + o(h^p)\) almost surely. This is the E&H power law with an amplitude \(a(x)\) that is a GP in \(x\).
  - It has the same structure as PRE 2401.07562 eq 10.
  - [lit] Bect §4.3: TWY2 intervals were "simultaneously smaller than the GCI interval and with a good coverage".
- **LB (lifted Brownian):** Boutelet & Sung 2503.23158 eq 4, with a variance envelope \(O(h^{2p})\).
  - [lit] Its parameter \(\gamma \in (0,1)\) "controls the correlation between increments". The Brownian kernel of Tuo–Wu–Yu, \(\min(h,h')^{2p}\), is the case \(\gamma = 0.5\), and only that case has independent increments.
  - [lit] For BM, Bect Prop. 2: the Richardson form "does not hold" a.s. Bect §2.2: the extrapolation "coincides with the observation of highest fidelity".
  - [lit] Boutelet Fig. 1: increments were positively correlated for an "average"-type QoI, and "somewhat uncorrelated or negatively correlated" for a "maximum"-type QoI. One fixed increment structure is therefore not generic, and LB learns it through \(\gamma\).
  - Our own evidence (Fig. B) is for BM, the case \(\gamma = 0.5\), only. LB with other values of \(\gamma\) is untested on our data.

Properties:
- The prior variance is \(Var[\delta(x,h)] = \sigma_\delta^2 k_x(x,x) h^{2p}\). It goes to 0 as \(h \to 0\) and is monotone in \(h\). This is the requirement "\(\sigma(x,h)\) decreases monotonically", applied to the discretisation error of a probe.
- The power law is **[assumption] A3**: one order \(p\) over \(X\) and an amplitude that varies with \(x\). It is NOT asserted for every QoI. Gate G0 tests it for each QoI and each transform. F1 and F2 show that it fails for \(\log \Delta t\) and is plausible for kLa.

### 2.4 Noise

$$ e_a \sim N(0, S_a), \quad S_a = D_a R_a D_a $$

- Runs are independent. \(D_a\) is the diagonal of the noise standard deviations \(s(x_k, h_a)\). \(R_a\) holds the correlations between the outputs of one run.
  - [lit] One run then counts as one correlated vector: Overstall & Woods 1506.04489 eq 2 (matrix normal). A separable form would force the noise correlation to equal the signal correlation, so \(R\) is kept separate.
- **Log-variance model.** \(\log s^2(x, h) = m_s + b_s \bar h + \zeta(x, h)\), where \(\zeta\) is a latent GP with variance \(\sigma_\zeta^2\) and Matérn-5/2 length scales \(\ell_\zeta\) in \(x\) and in \(\bar h\).
  - [lit] hetGP, 1611.05902 eq 14, smooths latent log-variances, so it works with few replicates.
  - The trend is \(b_s \bar h\), with a symmetric prior \(b_s \sim N(0, \tau_b^2)\) and \(\tau_b = \log 4\) [assumption]. One prior sd lets \(s\) change by a factor 2 between \(h_c\) and \(h = 0\). \(\bar h\) uses the fixed \(h_c\), so the prior is the same for every candidate structure. The spread may grow or shrink as \(h \to 0\).
  - The prior sd of \(m_s\) is 1.4 [assumption]: one sd is a factor 2 in \(s\).
  - [lit] Analogous evidence: in an LES study, the Lyapunov growth rate *increases* on finer meshes (1801.03046 §4.2: 285, 517, 589 per second). That is a turbulent flow, not ours, so it only motivates the symmetric prior.
  - [lit] Stroh 1605.02561 eq 3 and 1709.06896 eq 7e correlate the log-variances across levels; we keep that as the prior centre.
- **The latent field \(\zeta\) is part of the MCMC state**, not a plug-in. At a new \((x, h)\), \(\zeta\) is drawn from its GP conditional on the current sample.
- **Replicates.** The solver is deterministic, so a replicate is a run with a small perturbation (here, a shifted tracer-release cycle). \(s\) is the sensitivity of the QoI to that perturbation. Whether it equals the spread of a physical experiment is **untested** [assumption A7].
- \(s_0(x) = s(x, h^\star)\) is reported only if replicates exist at 3 or more levels, including the finest. Otherwise it is reported as "not identified".

### 2.5 Posterior sampling

- The state is \(\vartheta\): the hyperparameters \((\beta, p, \sigma_\delta, \ell_x, \ell_h\) or \(\gamma, \sigma_\mu, \ell_\mu, m_s, b_s, \sigma_\zeta, \ell_\zeta, \ell_\rho)\), plus the latent fields \(\zeta\) (about one value per distinct \((x, h)\), so about 400) and \(\tilde\mu\) if a link is used, plus the censored values.
- **Sampler: Gibbs with three blocks.**
  1. Latent Gaussian fields (\(\zeta\), \(\tilde\mu\)) are updated by **elliptical slice sampling**. [lit] Murray, Adams & MacKay 1001.0175 (abstract): it is for "models with multivariate Gaussian priors", "has no free parameters", and "works well for a variety of Gaussian process based models".
  2. The covariance hyperparameters of the latent fields are updated by the surrogate-data slice sampler of Murray & Adams. [lit] 1006.0868 (abstract): standard hyperparameter updates with non-Gaussian observations "require careful tuning and may converge slowly"; theirs "requires little tuning while mixing well in both strong- and weak-data regimes". The other hyperparameters use slice steps.
  3. Censored values are drawn from their truncated conditional (data augmentation) [inference].
  - HMC/NUTS is an alternative that I have not evaluated.
- **Convergence rule.** Run 4 chains. Use the draws only if the rank-normalised split-\(\hat R\) is below 1.01 and the effective sample size is above 400 for every reported quantity. [lit] Vehtari et al. 1903.08008 §1 ("only using the sample if \(\hat R < 1.01\)") and §4 (ESS "less than 400" indicates convergence problems). The Monte Carlo error of \(\sigma_{epi}\), a quantile spread, uses their quantile MCSE, and it must be below 5% of \(\varepsilon\) [assumption]. If the rule fails, run longer. Never report draws that fail it.
- Why sample at all: [lit] ML-II "underestimate[s] prediction uncertainty" (Lalchand & Rasmussen 1912.13440 §1), and MAP gave zero-width intervals in Stroh 1709.06896 §4.

### 2.6 Discrete structure: G0 filters, then stacking on an extrapolation score

- The structure index is \(M = (k_h\) family, \(g\), set of levels used\()\).
- **Order:**
  1. Gate G0 runs on every candidate \(M\). The candidates that fail it are removed. If none survive, the method stops (Section 4).
  2. **Extrapolation fits.** Fit every surviving \(M\) without the finest level, \(L_{max}\). The score is the density of the held-out \(L_{max}\) runs, so it tests extrapolation by one level in \(h\). This is where the families differ. Leave-one-run-out would mostly test interpolation, because the same \(u\) stays in at other levels.
  3. **Cross-fitting.** Split the \(L_{max}\) sites into two halves A and B, balanced in \(u\).
     - On A, compute stacking weights. [lit] Yao, Vehtari, Simpson & Gelman 1704.02030, eq 2.2 with the log score: \(\max_w \sum_i \log \sum_M w_M p(y_i \mid D_{-max}, M)\) with \(w \ge 0\) and \(\sum w = 1\).
     - **Regularisation.** We add the log of a Dirichlet(2, ..., 2) density to the objective, which pulls the weights toward equal. [lit] Yao et al. §4.1: with many models, "it can make sense to assign a strong prior … to the weights in estimation equation (2.2) to improve the regularization". The weight uncertainty is reported by a Bayesian bootstrap over the held-out sites [inference; Yao et al. use the Bayesian bootstrap for pseudo-BMA+].
     - **Low power.** With about 10 \(L_{max}\) sites, each half has about 5 sites against up to 8 structures. The weights are noisy, and G1 has little power. The report states both.
     - With those weights, run the calibration Gate G1 on B. Then swap A and B. G1 passes only if both folds pass. So the weights and the calibration check never use the same runs.
  4. **Final weights** for the prediction: stacking on all \(L_{max}\) runs. Then refit every \(M\) on all the data.
- **Why stacking.** [lit] Yao et al. (abstract): BMA "is flawed in the M-open setting", and they "recommend stacking of predictive distributions". Our G1 failures (Section 8) suggest that no candidate is true.
- **Why not selection.** In the level hold-out (Fig. B), the families differ by only about 2 nats in total log density, and their answers at \(h = 0\) differ 2–4×. Selecting one would hide that spread. (The 2 nats compare TWY2 with BM. LB with \(\gamma \neq 0.5\) is untested and can differ.)
- **Assumption A13:** weights that are good for extrapolating one level beyond the data are also good for extrapolating to \(h^\star\). This cannot be tested without a finer level.
- With 4 levels, the extrapolation fits use only 3 levels. This is the price of a held-out level, and it is a reason to add levels.
- **If only 2 levels remain** for the extrapolation fits (for example after G4 removes L6), \(p\) is not identifiable there. Then stacking is skipped, the weights stay equal over the G0 survivors (the prior), the output is flagged "weights prior-driven", and G1 cannot run, so the output is "uncalibrated" until a level is added.

### 2.7 Prediction on the physical scale (quantile based)

- For each structure \(M\) and each MCMC sample:
  - draw \(\mu(x) + \delta(x, h^\star)\) from the GP posterior, with \(\delta = 0\) at \(h^\star = 0\);
  - back-transform with \(g^{-1}\);
  - pool the draws over \(M\) with the weights \(w_M\).
- **Centre:** \(m_y(x)\) is the median of the pooled draws.
- **Spread:** \(\sigma_y(x) = (q_{84} - q_{16})/2\), from the 15.9% and 84.1% quantiles. For a Gaussian this equals the standard deviation.
- **Why quantiles.** Some back-transforms have no finite moments. If the rate \(z = -1/\Delta t\) is Gaussian, \(\Delta t\) has no finite mean or variance, because the density of \(z\) at 0 is not zero. Quantiles always exist.
- The Gaussian requirement is met by \(N(m_y, \sigma_y^2)\). Gate G6 compares the 2.5% and 97.5% quantiles with \(m_y \pm 1.96 \sigma_y\).
- With \(g = \log\), the posterior is lognormal before this match (question Q1).

### 2.8 Priors of every hyperparameter

All the priors are proper. The scales marked "§8" are problem-specific (R7).

| Parameter | Prior | Source |
|---|---|---|
| \(\beta\), no link | flat | universal kriging; the posterior is proper (Section 2.2) |
| \(\beta\), with link | \(N(b_{phys}, \mathrm{diag}(s_b^2))\); \(b_{phys}\) and \(s_b\) from §8 | required for a proper posterior |
| \((\sigma_\mu, \ell_\mu)\) | PC prior per input dimension: \(P(\ell < 0.1) = 0.05\), \(P(\sigma_\mu > S_\mu) = 0.05\); \(S_\mu\) from §8 | [lit] Fuglstad et al. 1503.00256 Thm 2.6 (joint PC prior for Matérn range and sd, \(d \le 3\)). Their range is \(\rho = \sqrt{8\nu}/\kappa\). With the Matérn form \(\kappa = \sqrt{2\nu}/\ell\), \(\rho = 2\ell\), so \(P(\ell < 0.1) = 0.05\) is \(P(\rho < 0.2) = 0.05\). [inference]: independent \(d = 1\) range priors per ARD dimension times one PC prior on \(\sigma\); and the use on the stationary factor of TWY2 |
| \((\sigma_\delta, \ell_x)\) | same PC form; \(P(\sigma_\delta > S_\delta) = 0.05\). \(\sigma_\delta\) is the error sd at \(h_c\), because \(\bar h = 1\) there. \(S_\delta\) from §8 | same |
| \(\ell_h\) (TWY2) | PC range prior, \(d = 1\): \(P(\ell_h < 0.1) = 0.05\) | same |
| \(\gamma\) (LB) | Uniform(0, 1) | [assumption] |
| \(p\) | \(\log p \sim N(0, 1)\) | G0 (Section 4) |
| \(m_s\), \(b_s\) | Section 2.4 and §8 | |
| \((\sigma_\zeta, \ell_\zeta)\) | PC prior: \(P(\sigma_\zeta > 1) = 0.05\), \(P(\ell_\zeta < 0.1) = 0.05\) | as above |
| \(\ell_\rho\) | PC range prior on the problem-specific output coordinate of \(R\) (§8), scaled to [0, 1] | as above |
| stacking weights | Dirichlet(2, ..., 2) penalty | Section 2.6 |

**Prior sensitivity (Gate G7).** Reweight (or refit) with each problem-specific scale (\(S_\mu\), \(S_\delta\), \(s_b\), the sd of \(\log p\)) halved and doubled, one at a time. The target is prior-dominated if, at any \(x \in \Sigma_N\), \(m_y\) moves by more than \(0.5 \sigma_{epi}\) or \(\sigma_{epi}\) changes by more than 20% [assumption]. A prior-dominated target is treated like a G0 failure: the method returns to the user.
- **Cost of G7:** 8 variants per surviving structure. Use importance reweighting of the base chains by the prior ratio [inference]. Accept the reweighted result only if the effective sample size of the weights is above 400. Otherwise run a full MCMC refit.

## 3. Bounded and reported quantities

All of these are on the physical scale, with \(\vartheta\) and \(M\) marginalised (Sections 2.6 and 2.7).

| Symbol | Definition | Type |
|---|---|---|
| \(\sigma_{epi}(x)\) | \(\sigma_y(x)\) of the target \(f(x, h^\star)\) | epistemic; more runs reduce it |
| \(s_0(x)\) | noise sd at \(h^\star\), back-transformed, with its band | aleatoric; compute cannot reduce it |
| \(\sigma(x,h)\) | sd of a new probe at \(h\): the posterior of \(f(x,h)\) plus the noise | the requirement's \(\sigma(x,h)\) |
| \(\sigma_{tot}(x)\) | the quantile spread of draws of \(g^{-1}(\mu + \delta + e)\) at \(h^\star\), with the noise added in \(g\)-space before the back-transform | the predictive for one new realisation; exact for any \(g\) |

- **The constraint bounds \(\sigma_{epi}\).** [lit] Design criteria use the de-noised variance because the noise is irreducible: Binois et al. 1710.03206 eq 2; 2412.07306 §2.1.
- A bound on \(\sigma_{tot}\) is infeasible if \(s_0 > \varepsilon\). The algorithm checks this (Step 4).
- **Deviation from requirement R5** (question Q1):
  - The *prior* discretisation variance is monotone in \(h\).
  - The *posterior* \(\sigma(x,h)\) is not monotone. [data] Fig. B: the bands are narrow at the levels that have data and wider at \(h = 0\). Any honest model must be less certain where no data exist. A monotone posterior would need a prior that forbids this.
  - \(s(x,h)\) has no imposed trend (Section 2.4).

---

## 4. Problem statement

**Given:**
- the prior \(\Pi\), the candidate structures \(M\), and the cost model (Section 6);
- \(\Sigma\), \(h^\star\), \(\varepsilon\), \(C\), and the batch size \(q\);
- the actions \(A\) = {new(\(u,h,T\)), extend(run \(j\), \(T \to T'\)), replicate(\(u,h\))}.

**P1.** Find a sequential policy that chooses batches \(a_1, ..., a_N\) and then stops, so that:

$$ \min \sum_i c(a_i) \quad \mathrm{s.t.} \quad \max_{x \in \Sigma_N} \sigma_{epi}(x \mid D_N) / \varepsilon(x) \le 1, \quad \sum_i c(a_i) \le C $$

and gates G0–G7 pass on \(D_N\). Here \(\varepsilon(x) = \varepsilon_{rel} m_y(x)\) (with \(m_y\) the median) or \(\varepsilon_{abs}\).

**P2.** If P1 is infeasible within \(C\): minimise \(\max_{\Sigma_N} \sigma_{epi}/\varepsilon\) subject to \(\sum c \le C\). Report "failed", the value reached, and where.
- A P2 result is reported as a prediction only if G0 and G1 pass. Otherwise it is reported as "uncalibrated".

Notes:
- [lit] Stacking designs (Sung, Ji, Mak, Wang, Tang, 2211.00268, eq 6 and 13) solves P1 for deterministic simulators with known cost. Ehara & Guillas 2104.02037 Prop. 2 solves the budget form P2.
- The bound is pointwise. It is not a simultaneous band.
- **Several QoIs from one run:** sum the criterion of Step 5 over QoIs \(j\) with weights \(1/\varepsilon_j^2\). [lit] Giles 1304.5472 §7.4 uses one Lagrange multiplier for several outputs.

### Gates. No probe at \(h = 0\) exists, so these are the only calibration checks.

| Gate | Test | Pass rule | Source |
|---|---|---|---|
| **G0 identifiability** (only if \(h^\star < h_{min}\); run on **each** candidate \(M\)) | Fit with the wide prior \(\log p \sim N(0, 1)\). Also report the observed increment ratios R. A known admissible range is in the prior (Section 2.2), so it is not tested here. | \(P(p > 0.5 \mid D) \ge 0.9\) [assumption] | [lit] E&H 2014 treat \(p < 0.5\) as anomalous; [data] Fig. D |
| G1 level hold-out | The cross-fitted test of Section 2.6, step 3: predict the held-out \(L_{max}\) runs of one half with the weights of the other half, including noise. Calibration only. | 95% coverage within binomial limits; sign test on z, p > 0.05; \(C_{LOO} = n^{-1}\sum z_i^2\) near 1 | [lit] Oliver 1311.0828 eq 17; Bachoc 1301.4320 eq 6 |
| G2 block LOO | Leave out one run (all its outputs), or one replicate set | z ~ N(0,1); Q–Q plot; the statistic \(U = \lvert I + E^T E \rvert^{-1}\) of Overstall & Woods eq 10 | [lit] Bachoc Prop. 3.1 (block form [inference]) |
| G3 noise | Replicate variance against the predicted \(s^2\) | chi-square test, p > 0.05 | [inference] |
| G4 pre-asymptotic | Refit without the coarsest level | the target moves less than \(\sigma_{epi}\) | [inference] |
| G5 structure | Posterior mean monotone along a coordinate that must be monotone (here \(\chi\)) | no violation on \(\Sigma_N\) | [lit] López-Lopera 1901.04827 eq 8 |
| G6 Gaussianity | Sample 2.5/97.5% quantiles against \(m_y \pm 1.96\sigma_y\) | difference < 10% of \(\sigma_y\) | [inference] |
| G7 prior sensitivity | Refit with each problem-specific prior scale halved and doubled (Section 2.8) | \(m_y\) moves < \(0.5\sigma_{epi}\) and \(\sigma_{epi}\) changes < 20% [assumption]; else treated like a G0 failure | [inference] |

**Repair rule.** When a gate fails, act in this order:
1. If G4 flags the coarsest level, remove it.
2. Recompute the stacking weights \(w_M\) (Section 2.6).
3. Buy a finer-level probe where \(|z|\) is largest.

Stop after **at most 2 repair cycles per gate**. If the gate still fails, the output is "uncalibrated". **If G0 fails for every candidate \(M\), the method stops and returns to the user** with five options (the same letters are used in Section 8 and Q3):
- (a) a finite \(h^\star\);
- (b) a finer level, with its cost forecast;
- (c) two fidelity indices;
- (d) a different QoI;
- (e) the E&H fallback. [lit] For \(p < 0.5\), E&H do not stop: they fit fixed orders (\(p = 1\), \(p = 2\), and the two-term form) with \(F_s = 3\). In our method, this is a prior concentrated on those orders, with its uncertainty reported. [data] For \(\Delta t\), per-rpm E&H gave \(U/\phi_0 \approx 0.9\)–\(1.4\) (diary 2026-10-02), which is far above \(\varepsilon = 10\%\).

**Deviation from R3.** E&H always return an estimate. Our default stops when G0 fails, and option (e) restores E&H behaviour on request.

\(\varepsilon\) is never widened without the user.

---

## 5. Relation to Yi et al. (2407.15110)

- **[lit] What Yi is.**
  - Eq 1: \(f^h(x) = g(f^l(x), x) + r(x)\). Eq 4 makes \(g\) polynomial.
  - The LF is a deterministic KRR, \(\rho\) comes from GLS, and \(r\) is a GP with homoscedastic noise (App. A eq 7).
  - It has two levels and no \(h\).
  - It is one step of NARGP (Cutajar 1903.07320 eq 2) and of RNA (Heo & Sung 2309.11772 eq 5).
- **Where it sits in our model.** The difference between two levels is

$$ g(y)(x, h_k) - g(y)(x, h_{k-1}) = \delta(x, h_k) - \delta(x, h_{k-1}) $$

  This is Yi's form with a linear \(g\), \(\rho_1 = 1\), \(\rho_0 = 0\), and a GP residual whose covariance follows from \(k_h\). A scale change \(\rho_1 \neq 1\) becomes an additive term under \(g = \log\).
- **What Yi cannot do:**
  - extrapolate to \(h^\star\) outside the data, because it has no \(h\);
  - carry the LF uncertainty, because the LF is a plug-in.
  - [inference] A chain of Yi steps with independent residuals has independent increments (the BM case), so Bect §2.2 applies.
- **Use.** Yi is a baseline in G1, and a fast approximation when one level has dense data.

---

## 6. Cost model

$$ c(a) = \kappa(u, h) \, (T_{spin}(a) + T_{obs}(a)) + c_0 $$

\(\kappa\) is the cost in core-s per simulated second:

$$ \log_2 \kappa = \kappa_0 + \gamma \ell + \gamma_2 (\ell - \bar\ell)^2 + \kappa_1(u) + \eta, \quad \eta \sim N(0, s_c^2) $$

- **Priors [assumption]:**
  - \(\kappa_0\): flat.
  - \(\gamma \sim N(3, 1)\), so the cost scales as \(2^{\gamma \ell}\), not \(2^\ell\) (question Q2). In 2D, a level gives 4× the cells and about 2× the time steps (`scripts/cost_model.py` measured 2.8×).
  - \(\gamma_2 \sim N(0, 0.5^2)\).
  - \(\kappa_1(u)\) is linear in problem-specific features of \(u\) with \(N(0, 1)\) coefficients (Section 8 uses \(\log rpm\)).
  - \(s_c\): half-normal with scale 0.5.
  - The fit is Bayesian linear regression.
- **[lit]** Guinet 2011.11456: simple low-variance cost models beat GPs. Snoek 1206.2944 §3.2 and taKG 1903.04703 §2.1 model the log cost.
- **[data]** \(\gamma\) is 2.9, 3.2 and 4.4 per level (Fig. C). That is why the quadratic term exists.
- \(c_0\) is the fixed overhead per run (start-up and I/O). Queue time is not counted.
- \(T_{spin}\) is the spin-up time. A warm start may reduce it, but **only after** test A11 passes (Section 8).
- **Cost that depends on the output.** [lit] Hutter 1310.1947 Def. 1: when the output is a time, the cost grows with the output. The expected observed time is

$$ E[T_{obs}] = \int_0^T P_n(\tau > t) \, dt $$

- extend(\(j\), \(T \to T'\)) costs \(\kappa (T' - T)\). [lit] taKG §2.6 prices a resumed run by the difference of costs.

![Fig. C](../experiments/multifidelity/cost_per_level.png)

*Fig. C. Measured cost, L6–L9, all rpm. Left: core-s per simulated second. Right: core-h to observe \(\Delta t_{0.95}\) after release, without the spin-up.*

**Does multilevel pay?**
- [lit] Stacking designs Thm 2 and Cor. 3 (eq 21 and 25): with bias rate \(\alpha\) and cost rate \(\beta\), multilevel beats single-level only if \(\alpha/\beta < 2\nu/d\). Most of the budget then goes to the fine levels (eq 23).
- For our case \(d = \dim(u) = 2\), because \(\chi\) comes free with each run. So \(2\nu/d = 2.5\) with \(\nu = 5/2\).
- The per-probe \(\beta\) is 4.1–5.6 (Fig. C, right).
- The check needs \(\alpha\) from data, so it is meaningful only after G0 passes. For mixing time at \(h^\star = 0\), it does not pass (Section 8).

---

## 7. Algorithm

**Step 0, priors (problem-specific).** Set the candidate \(M\), \(\phi\), the priors (\(p\), noise, cost), \(\Sigma\), \(h^\star\), \(\varepsilon\), \(C\) and \(q\).

**Step 1, initial design.**
- Use nested Sobol designs in \(u\), one per level. [lit] Stacking designs; its eq 9 gives the starting sizes.
- Necessary minimum [inference]:
  - at least 3 levels at 3 or more sites;
  - replicates at 3 or more levels at 3 or more sites.
- This minimum is not sufficient. [data] 4 levels × 10 rpm did not identify \(p\) for \(\log \Delta t\). G0 decides.

**Step 2, fit.** Run MCMC over \(\vartheta\) for each candidate \(M\).

**Step 3, gates.** Run G0–G7, then apply the repair rule.

**Step 4, stop and forecast.**
- Success: P1 holds and the gates pass.
- If the bound is on \(\sigma_{tot}\) and \(\max s_0 > \varepsilon\): infeasible. Stop and report.
- **Forecast.** Run the greedy policy of Step 5 on fantasy data until success. Repeat over posterior samples. This gives a distribution of the cost to success. If the probability that spent plus forecast cost exceeds \(C\) is above 0.5, report "probably infeasible" before buying more runs.

**Step 5, acquisition: maximum rate of uncertainty reduction per cost.**

$$ a^\star = \arg\max_{a \in A} \; (H_n - E J_n(a)) / E c(a), \qquad H_n = \sum_{x \in \Sigma_N} (\sigma_{epi}^2(x) / \varepsilon^2(x) - 1)_+ $$

- [lit] The ratio is MR-SUR (Stroh 2007.13553 eq 11). In that paper, \(H_n\) is the integrated posterior variance (eq 6).
- The hinge form above is **my change** [inference]. It is zero exactly when P1 holds on \(\Sigma_N\). It is not smooth, but \(a\) is chosen from a finite candidate set, so this does not matter.
- **\(J_n(a)\) is a Monte Carlo average over fantasies, not a deterministic value.** Per fantasy:
  - (i) Draw \(\vartheta\) from the posterior, including the noise field \(\zeta\).
  - (ii) Draw the outputs of \(a\) from the predictive. This also draws which nested outputs are reached before \(T\) and which are censored.
  - (iii) Update the posterior, reweighting the \(\vartheta\) samples by the fantasy likelihood, and compute \(H\).
  - [lit] Snoek 1206.2944 §3.3 uses Monte Carlo over pending outcomes.
  - Reweighting gives value to probes that reduce the uncertainty of \(p\).
- **Run length.** \(T\) is the 0.95 predictive quantile of the time to the last needed output. [lit] taKG §2.3 charges the cost at the largest trace point; the shorter outputs come free.
- **Batch of \(q\).** Choose greedily, and condition each choice on the pending runs. [lit] Takeno 1901.08275 eq 7–8; Kathuria 1611.04088.

**Computational cost of Steps 2–5 [inference, order of magnitude].**
- Mixing example: about 40 runs × 10 outputs = 400 data. One GP solve costs about \(400^3/3 \approx 2 \times 10^7\) flops.
- MCMC: \(10^4\) samples × 4 structures, which is minutes on one node.
- Acquisition: about 120 candidates (10 \(u\) × 4 levels × 3 action types) × 16 fantasies × 200 thinned samples. Each needs a rank-10 update of \(O(n^2 m) \approx 2 \times 10^6\) flops. That is about \(8 \times 10^{11}\) flops, which is minutes to an hour on one node. This is negligible against one L8 run.
- **Quantile spreads:** for each candidate and fantasy, about 800 pooled draws (4 structures × 200 samples) at about 300 points of \(\Sigma_N\), then a sort: about \(10^7\) operations, so \(2 \times 10^{10}\) in total. This is small against the GP updates.
- **Number of fantasies \(F\)** [assumption]: start with \(F = 16\). Double \(F\), up to 256, until the Monte Carlo standard error of \(H_n - E J_n(a)\) for the best candidate is below 10% of the gap to the second-best candidate. The quantile spread is robust to heavy tails, which helps.
- **MCMC** (Section 2.5): latent fields make each chain slower. The time must be measured, and the estimate above (minutes) is only for the hyperparameters.
- **G7:** 8 reweightings per surviving structure, which is cheap; a full refit only where the weight ESS is below 400.
- **Forecast (Step 4):** a greedy plan of about 20 steps at the full cost would be about 20 acquisitions. To keep it cheap, the forecast uses 20 fixed \(\vartheta\) samples, no reweighting, and the expected reached set. This is an approximation, and the realised cost is compared with it after each batch.

**Step 6, run.**
- If an output is not reached by \(T\), extend the run from its checkpoint. [lit] Freeze-thaw, Swersky 1406.3896 Alg. 1.
- Otherwise keep the censored value \(\tau > T\) in the likelihood. [lit] Hutter 1310.1947: a capped run gives a lower bound. Data augmentation in the MCMC is [inference].
- Go to Step 2.

**Output.** For each \(x \in \Sigma\): \(m_y\), \(\sigma_{epi}\), \(s_0\) (or "not identified"), and \(\sigma_{tot}\). Also: the stacking weights \(w_M\), the posterior of \(p\), the gate table, the cost and the allocation per level.

**Novelty check.** [lit, negative] No paper read in this session combines all four of these: a Richardson-type or LB kernel with an \(x\)-dependent amplitude; a marginalised \(p\); heteroscedastic noise that is correlated inside a run; and nested outputs whose cost depends on the output.

---

## 8. Instance: mixing time \(\Delta t_{mix}((\chi, \theta_{max}, rpm), h)\)

| Item | Choice | Status |
|---|---|---|
| \(x\), \(u\) | \(x = (\chi, \theta_{max}, rpm)\), with \(\chi\) in 0.50, 0.55, ..., 0.95; \(u = (\theta_{max}, rpm)\) | |
| \(\Sigma\) | χ in [0.5, 0.95], rpm in [15, 37.5]; the \(\theta_{max}\) range is open (Q2) | |
| \(g\) candidates | \(\log\), \(-1/\Delta t\) | [data] for predicting the next level, log wins by 6–16 nats (Fig. B); for convergence, only the rate is monotone (Fig. D) |
| \(\varepsilon\) | \(\varepsilon_{rel} = 0.10\) | user |
| \(\phi\) | \((1, \log(-\ln(1-\chi)), \log rpm, \log \theta_{max})\) | A single exponential gives the second term. Fig. A shows that \(\lambda\) changes 0.3–3.7× along a curve, so this is a mean trend only. |
| \(p\) | the G0 prior \(\log p \sim N(0, 1)\) | [assumption] |
| noise | \(m_s \sim N(\log 0.03^2, 1.4^2)\), \(\tau_b = \log 4\) | [assumption]; replicates must replace it |
| \(S_\mu\), \(g = \log\) | 1 (\(\mu\) varies by up to a factor \(e^2\) over \(\Sigma\) at 2 sd) | [assumption] |
| \(S_\delta\), \(g = \log\) | 4 (the error at L6 can be a factor 50) | [assumption], **set after seeing L6–L9** (×40–60 from L6 to L9) and justified by the Péclet argument of F1; G7 tests its influence |
| rate link: \(b_{phys}\), \(s_b\) | \(\tilde\mu = \log\) rate. Intercept: \(\Delta t_{0.95}\) = 30 rocking periods at 25 rpm, with \(s_b = 2\) (a factor 7.4 per sd). **No source:** this is a guess, and it conflicts with the data (L9 at 25 rpm is 304 s, about 126 periods, and every ladder grows with refinement). So it leans toward a short, finite limit. G7 tests its influence. A sourced value is needed (Q1). Slopes: −1 on \(\log(-\ln(1-\chi))\) (single exponential), +1 on \(\log rpm\) (time in periods), 0 on \(\log \theta_{max}\); each with \(s_b = 1\). | [assumption]; Kim's values are **not** used, so the comparison with Kim stays independent |
| \(S_\mu\), \(S_\delta\), rate link | \(S_\mu = 1\) on \(\tilde\mu\); \(S_\delta = 40 r_{ref}\), with \(r_{ref}\) the rate of the \(b_{phys}\) intercept | [assumption], set after seeing the data, as above |
| \(R\) coordinate | \(\psi = \log(-\ln(1-\chi))\), scaled to [0, 1] over \(\Sigma\) | |
| cost features | \(\kappa_1(u) = k \log rpm\) | [data] cost per cycle changes 4× across rpm at L10 (`scripts/cost_model.py`) |
| \(R\) | \(\exp(-\lvert \psi_k - \psi_{k'} \rvert / \ell_\rho)\) | [assumption]; G2 tests it |
| \(\sigma^2_{max}\) | the tracer variance at release; exact value 0.25 (top-half release); a run is accepted if it is within ±20%, that is [0.20, 0.30] | `scripts/postprocess.py` (`_SIGMA2_MAX_TOL = 0.20`) |

**Gate status with the data we have (per-rpm fits, no pooling; `scripts/diag_h_kernel_test.py`, `scripts/diag_observed_order.py`):**
- **G0 at \(h^\star = 0\)** (provisional: per-rpm triplets and the parametric pooled rate fit; G0 under the GP kernels is **untested** and is action 1):
  - \(g\) = identity: R = 0.2–0.4, so the increments grow. Expected to fail.
  - \(g = \log\): R ≈ 1 (0.6–1.2), so there is no convergence. Expected to fail.
  - Level sets: removing L6 does not help. On L7–L9 alone, \(\log\) has median R = 0.89–1.24.
  - \(g\) = rate, with positivity in the prior (\(\mu = -\exp(\tilde\mu)\)): the parametric pooled fit identifies \(p = 1.41 \pm 0.31\), so it may pass G0.
    - But the data put \(r_\infty\) near 0: the per-rpm Richardson limits are negative at 9 of 22 triplets. With positivity in the prior, the posterior of the rate then piles near 0, and the upper quantiles of \(\Delta t\) become very large.
    - So the method goes on to the Step 4 forecast. That forecast is expected to report "probably infeasible" for \(\varepsilon_{rel} = 0.10\) [inference; action 1 computes it].
  - The verdicts depend on two choices: positivity in the prior, and the G0 threshold.
- **G1 fails.** In log space, 29/29 BM z-scores and 28/29 TWY2 z-scores are positive (fit L6–L8, predict L9). This is consistent with F1.
- **[data]** L6 loses 7–10% of the tracer mass (L7 3%, L8 0.7%, L9 under 0.4%). So G4 will probably remove L6.

![Fig. B](../experiments/multifidelity/h_kernel_test_chi95.png)

*Fig. B. \(\Delta t_{0.95}\) against cell size, per rpm. TWY2 and BM posteriors fitted on L6–L9, with 95% bands. The cross is the single L10 run. The bands between levels are not used.*

**The L10 anomaly.**
- fig9_l10b_seg2 at 32.5 rpm gives 18.7 / 31.6 / 64.7 s, against 39.5 / 72.9 / 128.3 s at L9.
- Its \(\sigma^2_{max} = 0.20\) is at the edge of the tolerance.
- It was a warm-start chain with 33 cycles before release. L9 was a cold start with 80 cycles.
- It is not used until test A11 removes the protocol confound.

**Options when G0 fails or the forecast says infeasible** (the user decides, Q3; the letters match Section 4):

| Option | What it is | Status |
|---|---|---|
| (a) Finite target \(h^\star\) | Inside or one level beyond the data, G0 is not needed or is weaker. | Fig. B: the BM band at L10 is still 56–290 s, so more L9/L10 data are needed even for this. Also: the upstream Kim code sets `MAXLEVEL = 9` (`tests/fixtures/kim_upstream/BioReactor.c:180`), but our notes call L10 "Kim's grid". This must be resolved before we choose \(h^\star\). |
| (b) A finer level | Add L10 runs. | [data] Per probe, the cost to \(\Delta t_{0.95}\) grew 50× from L8 to L9 (Fig. C); one more level is of the order of \(10^4\) core-h per probe [inference, extrapolated]. |
| (c) Two fidelity indices | Use \(h = (h_{flow}, h_{scalar})\): the passive tracer on a finer grid, on a replayed flow period. | [lit] CONFIG eq 19 has a kernel for several fidelity parameters. The replay solver is in BACKLOG and is untested. |
| (d) A different QoI | A mixing measure that is less sensitive to diffusion. | Not read; I would read sources first. |
| (e) E&H fallback | Fixed orders with \(F_s = 3\) (Section 4). | [data] per-rpm \(U/\phi_0 \approx 0.9\)–\(1.4\), far above 10%. |

**First actions, in order:**
1. **Pooled fit, no compute.** Run G0 and G1 on all L6–L9 × 10 rpm × 10 χ values (plus the `kmix_l7_th2`–`th6` angle runs, so that \(\theta_{max}\) varies), with the generic model, for \(\Delta t\) (log, rate) and kLa (identity, log). Pooling may identify \(p\) where per-rpm fits cannot. If G0 still fails, F1 is confirmed for the generic model.
2. **Test A11, about 20 core-h.** L8 at 32.5 rpm with a 33-cycle spin-up, against the existing 80-cycle run. This decides whether the L10 point is usable and whether warm starts are allowed.
3. **Replicates, about 100 core-h at L7/L8.** Add about 600 core-h for one L9 pair. Shift the release by a few cycles. This gives \(s(x,h)\) at 3 levels, its trend, and \(R\).

---

## 9. Assumption register

| # | Assumption | Status | Test |
|---|---|---|---|
| A1 | \(g(y)\) is Gaussian | untested | G6; Q–Q plots of the G2 residuals |
| A2 | The code converges, \(\delta(x,0) = 0\), within the affordable levels | [data] **fails** for \(\Delta t\) on L6–L9 (F1); plausible for kLa (F2) | G0 |
| A3 | Leading power law with one \(p\) over \(X\) | [data] it fails for \(\log \Delta t\); [lit] Bect Prop. 3 under TWY2 | G0; then let \(p\) vary with \(x\) |
| A4 | \(k_h\) family | [data] TWY2 vs BM (LB with \(\gamma = 0.5\)) not separable yet (~2 nats); LB with other \(\gamma\) untested | G2 weights, G1 |
| A5 | Separable \(k_x k_h\) | untested | G2 residuals against \(x\) |
| A6 | Noise trend in \(h\) of either sign | [lit] analogy only | replicates (action 3) |
| A7 | A release-shift replicate represents physical spread | untested | experimental replicates (outside scope) |
| A8 | Stationary kernels in rpm | untested | G2 residuals near 22.5 rpm |
| A9 | Log-cost quadratic in level | [data] linear fails (γ = 2.9/3.2/4.4) | posterior predictive check on new runs |
| A10 | The coarsest level is asymptotic | [data] doubtful for L6 | G4 |
| A11 | Spin-up length does not change \(\Delta t\) | untested; the L10 point conflicts with it | action 2 |
| A12 | Fantasy reweighting values \(p\)-reducing probes correctly | [inference] | synthetic test |
| A13 | Stacking weights for one-level extrapolation transfer to \(h^\star\) | untestable without a finer level | a finer level |

---

## 10. Questions for the user

- **Q1.**
  - May the constraint bound \(\sigma_{epi}\) and report \(s_0\) separately?
  - Do you accept that the posterior \(\sigma(x,h)\) is not monotone (Fig. B), while the prior discretisation variance is?
  - Do you accept the Gaussian \(N(m_y, \sigma_y^2)\) built from the median and the 68% quantile spread, instead of the mean and the variance (Section 2.7)? With \(g = \log\), the underlying posterior is lognormal.
  - Do you accept averaging over structures (Section 2.6) instead of choosing one?
  - Some prior scales (\(S_\delta\), and \(S_\mu\) for the rate) were set after seeing L6–L9. That uses the data twice, a small deviation from fully Bayesian; G7 tests its influence. Do you accept this, or can you give physical values? The same applies to the rate-link intercept (30 periods), which has no source.
- **Q2.** What are the \(\theta_{max}\) range, the budget \(C\), and the batch size \(q\)? You wrote "cost scales with \(2^N\)". The data say \(2^{\gamma N}\) with \(\gamma\) = 2.9–4.4 per level (Fig. C). May \(\gamma\) be learned?
- **Q3.** For mixing time at \(h = 0\), G0 fails or the forecast is expected to say infeasible (Section 8). Which option do you want: (a) a finite \(h^\star\), (b) a finer level, (c) two fidelity indices, (d) a different QoI, or (e) the E&H fallback? And which grid did Kim use: L9 (the upstream code) or L10 (our notes)?
- **Q4.** Do you accept Yi et al. as a baseline and special case (Section 5)?
- **Q5.** Is model-form error (2D vs 3D) in scope? If yes, which experimental data can we use?
- **Q6.** May I run actions 2 and 3 (about 120 core-h, or about 720 with the L9 pair)?
