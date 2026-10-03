# Yi-h: grid-convergent multi-fidelity learning (problem definition, algorithm, assumptions)

Draft 12, 2026-10-03. This draft makes the method an extension of Yi et al. (2024), generic for engineering QoIs, with the bioreactor as a worked example only.

**History.** Drafts 2–6 went through 5 rounds of adversarial review; the reviewer approved draft 6. Draft 8 rewrote the structure after the user's direction of 2026-10-03. A fresh reviewer rejected draft 8 (6 major holes) and draft 9 (2 major holes), and approved draft 10 with 6 minor holes. Draft 11 closed those. Draft 12 makes the model fully Bayesian (user direction), rewrites the literature position, and adds the spin-up test.

**Status labels:**
- **[lit]**: I read it in the cited paper.
- **[data]**: our measurement; the script is named.
- **[inference]**: my reasoning; it is not tested.
- **[assumption]**: a choice; it is not tested.

---

## 0. Summary

**The engineering problem.** A simulation code computes a quantity of interest (QoI) \(y(x, h)\) at inputs \(x\) and numerical resolution \(h\). Examples are mean or maximum wall shear, a mass-transfer coefficient, or a mixing time.
- Finer resolution costs much more (often \(2^{\gamma}\) per halving of \(h\), with \(\gamma\) from 2 to 5).
- The engineer wants the converged value \(f(x)\) over a design region, with a stated uncertainty, inside a compute-time budget. \(f(x)\) is the median outcome of a run as \(h \to 0\). This equals the mean when the run-to-run spread is symmetric on the physical scale (Section 2.2 gives the reason for the median).
- Two literatures each solve half of this problem (details and sources in the guide, Section 11):
  - **Solution verification** (Eça & Hoekstra 2014; probabilistic Richardson extrapolation, 2401.07562) estimates the converged value at **one** condition.
  - **Multi-fidelity regression** (Yi et al. 2407.15110, and the others read) predicts over **all** conditions, but its target is the most expensive fidelity that was run, not \(h = 0\).
  - A small line of work joins the two halves (Tuo–Wu–Yu; CONFIG 2209.13748; Boutelet & Sung 2503.23158; DNA 2506.08328; Stroh et al.). In the papers read, none combines a learned rate, calibrated uncertainty at \(h = 0\), and sequential budget-aware design over a design region.

**The method (Yi-h).** It keeps the structure of Yi et al.'s KRR-LR-GPR: abundant cheap data, a linear transfer between fidelities, and a Bayesian residual. It makes every part Bayesian, and it adds four things:
1. Any number of resolution levels in one likelihood. No level is the truth.
2. A covariance in \(h\) that encodes convergence (an E&H-type power law, \(h^p\), with \(p\) learned).
3. The target is \(h = 0\). The uncertainty of that extrapolation is part of the answer.
4. Cost-aware sequential design. It decides which inputs and which resolutions to run next, until the uncertainty is below a tolerance or the time budget is used.

Optional modules use special structure when a problem has it (Section 4). The bioreactor (Section 6) uses some of them.

---

## 1. Problem definition (generic)

**Given:**
- **Inputs:** \(x \in X \subset R^d\), with \(x = (u, v)\). Here \(u\) are the controls that a run sets (for example rpm and angle), and \(v\) are output coordinates that one run returns many of (for example a threshold χ). \(v\) can be empty.
  - \(\Sigma \subset X\) is the region of interest. \(\Sigma_N\) is a scrambled Sobol set of \(100 d\) points in \(\Sigma\), on which every "max over \(\Sigma\)" is computed.
  - All coordinates are scaled to [0, 1] over \(\Sigma\) inside the kernels.
- **Resolution:** a vector \(h = (h_1, ..., h_k) \in (0, \infty)^k\), with \(k \ge 1\). Examples: a cell size, a time step, the cell size of a scalar solved on its own grid.
  - Each component takes discrete values with any refinement ratio, \(h_{j,\ell+1} = h_{j,\ell}/r_{j,\ell}\) with \(r_{j,\ell} > 1\). Finer values are always allowed; the budget is the only limit.
  - Inside the kernels, \(\bar h_j = h_j / h_{j,c}\), where \(h_{j,c}\) is fixed at Step 0 (the coarsest candidate value). It does not change when a level is removed later, so the priors keep their meaning.
  - Vectors are ordered componentwise: \(h \le h'\) when \(h_j \le h_j'\) for every \(j\).
  - \(h = 0\) means the exact solution of the model equations.
- **Probe:** one run \(a = (u, h, T)\), with run length or other run setting \(T\). It returns a vector \(y_a\) at the coordinates \(O(a) \subset X\). In the generic case \(O(a) = \{u\}\).
- **Cost:** \(c(a) > 0\), in core-hours. It is unknown and is learned.
- **Budget** \(C\) (core-hours) and **tolerance** \(\varepsilon\) on the physical scale, relative (\(\varepsilon_{rel}\)) or absolute.

**Quantities to model and report** (the uncertainties are defined in Section 2.8):

| Symbol | Meaning | Decreases as \(h \to 0\)? |
|---|---|---|
| \(f(x)\), \(m_y(x)\) | the converged value (median run outcome at \(h = 0\)), and its posterior median | — |
| \(\sigma_{epi}(x)\) | epistemic: uncertainty about \(f(x)\); more runs reduce it | it decreases as data are added |
| \(s_0(x)\) | aleatoric: run-to-run spread at \(h = 0\) | no trend is imposed (Section 2.4) |
| \(\sigma_{env}(x, h)\) | fidelity envelope, in \(\Lambda\) units: the RMS size of the error of level \(h\) under the posterior error model, mean and random part (Section 2.8) | **yes**, monotone to 0 by construction (Fig. E, curve 1) |
| \(\sigma_{fid}(x, h)\) | posterior RMS error of level \(h\) against \(f(x)\), on the physical scale | it is 0 at \(h = 0\) and tends to the true error \(\lvert\delta(x,h)\rvert\), so it is monotone only if the true convergence is monotone (Fig. E, curve 2) |
| \(\sigma_{know}(x, h)\) | knowledge: the posterior sd of \(f(x, h)\) | **no**: it is small where data exist (Fig. E, curve 3) |

**Problem P1.** Find a sequential policy that chooses batches of probes \(a_1, ..., a_N\) and stops, so that

$$ \min \sum_i c(a_i) \quad \mathrm{s.t.} \quad \max_{x \in \Sigma_N} \sigma_{epi}(x \mid D_N) / \varepsilon(x) \le 1, \quad \sum_i c(a_i) \le C $$

and the calibration gates pass (Section 3). Here \(\varepsilon(x) = \varepsilon_{rel} m_y(x)\) or \(\varepsilon_{abs}\).

**Problem P2,** when P1 is infeasible within \(C\): minimise \(\max_{\Sigma_N} \sigma_{epi}/\varepsilon\) subject to \(\sum c \le C\). Report "not met", the value reached and where.

Notes:
- The constraint bounds \(\sigma_{epi}\), because compute can reduce only that part. \(s_0\) and \(\sigma_{tot}\) are modelled and reported (both defined on the physical scale in Section 2.8).
  - [lit] Design criteria use the de-noised variance for this reason: Binois et al. 1710.03206 eq 2.
- The bound is pointwise. It is not a simultaneous band.
- [lit] Stacking designs (Sung, Ji, Mak, Wang, Tang, 2211.00268 eq 6, 13) solves P1 for deterministic codes with known cost. Ehara & Guillas 2104.02037 Prop. 2 solves the budget form P2.

![Fig. E](../experiments/multifidelity/sigma_concepts.png)

*Fig. E (synthetic, `scripts/plot_sigma_concepts.py`). Left: data at 4 levels and the posterior of \(f(h)\). Right, against \(h\):*
- *(1) the fidelity envelope \(\sigma_{env}\) (here \(\Lambda\) is the identity, so its units equal the physical scale);*
- *(2) the posterior RMS error of level \(h\), which is \(\sigma_{fid}\);*
- *(3) the posterior sd of \(f(h)\), which is \(\sigma_{know}\);*
- *(4) the sd of \(f(0)\) after one more probe at \(h\).*

*(1) goes to 0 monotonically by construction: "finer is closer to the truth" is in the prior. (2) is monotone here because the synthetic truth converges monotonically. With an oscillating truth, \(1 + 0.4 h^{1.5} \cos(2\pi h)\), the reviewer measured (2) = 0.186 at \(\bar h = 0.5\) and 0.117 at 0.75. Oscillatory convergence is real (Fig. D, 37.5 rpm). (4) decreases for finer probes: "finer is more informative". (3) is not monotone, because it describes what we know about each level, and we know most where we measured.*

---

## 2. Model: Yi-h

### 2.1 What Yi et al. do [lit, 2407.15110]

- Eq 1: \(f^h(x) = g(f^l(x), x) + r(x)\). The three parts are an LF model trained on LF data, a transfer model, and a residual model.
- Eq 2 is the linear special case \(f^l(x)\rho + r(x)\). Eq 4 allows a polynomial \(g\).
- The LF model is a deterministic KRR. The coefficients \(\rho\) come from GLS (Algorithm 1). \(r\) is a GP, with homoscedastic noise (App. A eq 7).
- The setting is "scarce, resource-intensive high-fidelity data with abundant but less accurate low-fidelity data" (abstract).
- There are two levels and no resolution variable. The target is \(f^h\).

### 2.2 The extension: Yi's transfer as a function of resolution

User direction (2026-10-03): be more Bayesian than Yi, and accept the cost; the problems have at most about 10 inputs. So no part of the model is a deterministic plug-in.

For output \(k\) of probe \(a\):

$$ \Lambda(y_{a,k}) = \rho_0(x_k, h_a) + \rho_1(x_k, h_a) \, \mu(x_k) + \delta(x_k, h_a) + e_{a,k} $$

$$ \rho_0 = \sum_{j=1}^{k} c_{0j} \bar h_j^{p_j(x)}, \qquad \rho_1 = 1 + \sum_{j=1}^{k} c_{1j} \bar h_j^{p_j(x)} $$

| Term | Role | Relation to Yi |
|---|---|---|
| \(\mu(x)\) | the converged value in \(\Lambda\)-space; a GP with a linear mean \(\beta_0 + \beta^T x\) (or a physical basis) | Yi's high-fidelity target, moved to \(h = 0\) |
| \(\rho_0, \rho_1\) | a linear transfer between the converged value and level \(h\). The coefficients go to \((0, 1)\) as \(h \to 0\) with the E&H power law. | Yi's linear transfer (eq 2, eq 4 with \(M = 2\)), made a function of \(h\) |
| \(\delta(x, h)\) | GP residual of level \(h\), with the convergence kernel of Section 2.3; \(\delta(x, 0) = 0\) | Yi's residual \(r\), made a function of \(h\) |
| \(e_a\) | run-to-run noise vector, heteroscedastic, correlated inside a run | generalises Yi's homoscedastic noise |
| \(\Lambda\) | monotone output transform (identity, log, reciprocal, ...). Named \(\Lambda\) because Yi's \(g\) is the transfer model. | new; chosen per QoI (Section 2.7) |

- **Every level is in one likelihood**, the cheapest included. \(\mu\), \(\rho\), \(\delta\) and all hyperparameters are inferred jointly (Section 2.5).
- **Exact inference:** for given \((c, p)\) and hyperparameters, the model is linear and Gaussian in \((\mu, \delta)\). So these two integrate out exactly, and the MCMC runs only over the low-dimensional rest.
- **Cost:** at most about 10 inputs and a few thousand runs, so exact GP algebra is affordable. One solve with \(n = 2000\) costs about \(3 \times 10^9\) flops.
- **Target:** \(f(x) = \Lambda^{-1}(\mu(x))\).
  - The noise is symmetric with zero median in \(\Lambda\)-space, and \(\Lambda\) is monotone. So \(f(x)\) is the **median** run outcome at \(h = 0\) for **every** transform, because a median does not change under a monotone map. The user accepted the median (2026-10-03).
  - So all candidate structures estimate the same quantity, and Section 2.7 can pool them.
- **Optional external cheap source** (module S8): a different, cheaper model, not a resolution of this code (a reduced-order model, a correlation).
  - Its data are \(y_s = f_s(x) + e_s\), with \(f_s\) a GP, and \(\mu(x) = \beta_s f_s(x) + r(x)\), all in the same likelihood.
  - [lit] This is Yi eq 2 between the source and the converged value, with a GP for the low fidelity. Yi §2 describes this as the data-scarce literature's choice: "MF models use GPRs for both \(f^l(x)\) and \(r(x)\) and implicitly assume a linear transfer-learning model". Yi's own choice for that slot is a deterministic KRR.

**Yi is a special case [inference, by construction].** Take two levels \(h_l > h_h\) with \(\Lambda\) the identity, and remove \(\mu\) from the two equations:

$$ f(x, h_h) = \rho_0' + \rho_1' \, f(x, h_l) + r'(x) $$

- Here \(\rho_1' = \rho_1(h_h)/\rho_1(h_l)\) and \(\rho_0' = \rho_0(h_h) - \rho_1' \rho_0(h_l)\).
- \(r' = \delta(x, h_h) - \rho_1' \delta(x, h_l)\) is a GP. It contains the coarse level's own error \(\delta(x, h_l)\), so it is **correlated with** \(f(x, h_l)\). Yi's residual is independent of the low-fidelity predictor. So the result has Yi's **form**, and it is Yi's **model** only if \(\delta(x, h_l) \equiv 0\), or if \(r'\) is uncorrelated with \(f(x, h_l)\), which needs \(k_h(h_h, h_l) = \rho_1' k_h(h_l, h_l)\) (a Markov property in \(h\) that does not hold in general). Otherwise Yi's fitted \(\rho\) absorbs part of the grid error (a fresh reviewer found this on the guide's toy: correlation \(-0.93\)).
- This is Yi's eq 2 (and eq 4 with \(M = 2\)) when that holds and all of these also hold:
  - the orders \(p_j\) do not vary with \(x\), so \(\rho'\) is constant;
  - \(f(x, h_l)\) is replaced by a deterministic KRR fit (Yi's LF model);
  - the residual kernel is Yi's RBF (Yi eq 5);
  - the noise is homoscedastic and uncorrelated (\(R = I\), \(\zeta \equiv 0\), \(b_s = 0\));
  - the coefficients have a flat prior;
  - the hyperparameters come from ML-II with the concentrated likelihood (Yi Algorithm 1, step 2.2);
  - the target is \(f(x, h_h)\), not \(f(x, 0)\).

**What Yi-h adds to Yi:**
- the transfer coefficients and the residual are functions of \(h\) that converge to the identity and to zero;
- any number of levels in one likelihood;
- the target at \(h = 0\);
- full Bayesian inference, including the order \(p\);
- heteroscedastic, within-run-correlated noise;
- cost-aware design.

### 2.3 Convergence covariance in \(h\)

For one resolution component:

$$ \delta \sim GP(0, \; \sigma_\delta^2 k_x(x, x') k_h(h, h')) $$

Two families are candidates for \(k_h\), and Section 2.7 weights them. Several components and an order that varies with \(x\) follow below.

- **TWY2 (Richardson type)**, for one component: \(k_h = (\bar h \bar h')^{p} c_\nu(\bar h - \bar h'; \ell_h)\), with \(c_\nu\) a Matérn correlation with \(\nu = 3/2\). [lit] Bect §4.3 found \(\nu = 1/2\) or \(3/2\) good. Our diagnostic scripts used 5/2.
  - [lit] Bect et al. 2103.14559 Prop. 3 (one component): \(\delta(h) = A h^p + o(h^p)\) almost surely. This is the E&H power law, with the amplitude a GP in \(x\).
  - [lit] It has the same structure as PRE 2401.07562 eq 10.
  - [lit] Bect §4.3: TWY2 intervals were "simultaneously smaller than the GCI interval and with a good coverage".
- **LB (lifted Brownian):** [lit] Boutelet & Sung 2503.23158 eq 4, with \(\gamma \in (0,1)\) controlling "the correlation between increments".
  - The Brownian kernel of Tuo–Wu–Yu, \(\min(h,h')^{2p}\), is \(\gamma = 0.5\).
  - [lit] Bect Prop. 2: in that case the Richardson form "does not hold".
  - [lit] Boutelet Fig. 1: increments are positive for an average-type QoI, and "somewhat uncorrelated or negatively correlated" for a maximum-type QoI.
- **Several resolution components: an additive error.**

$$ \delta(x, h) = \sum_{j=1}^{k} \delta_j(x, h_j), \quad \delta_j \sim GP(0, \; \sigma_{\delta,j}^2 k_{x,j}(x, x') k_{h,j}(h_j, h_j')) $$

  - The \(\delta_j\) are independent. The error vanishes only when **every** component goes to 0. A product kernel would vanish when any one component goes to 0, so refining the flow grid alone would remove the scalar-grid error, which is wrong.
  - [lit] Boutelet & Sung §2.1, citing Ji et al.: the error must stay non-negligible while any component is nonzero. CONFIG 2209.13748 eq 19 is an alternative.
- **Order that can vary with \(x\)** (assumption A3; **TWY2 only**): \(\delta_j = \bar h_j^{p_j(x)} e_j(x, h_j)\) with \(e_j\) a stationary GP, and \(\log p_j(x) = \log p_{j0} + \pi_j(x)\). LB keeps one shared order per component, because its order sits inside the kernel's power and no varying form was read for it.
  - \(\pi_j\) is a GP with a PC prior that shrinks its variance to 0, so the base model is one shared order [inference]. The data switch the variation on only if they need it.
  - The covariance stays valid: \(b(z) b(z') k(z, z')\) is positive semi-definite for any function \(b\), here \(b = \bar h^{p(x)}\).
- **Properties:**
  - The prior variance of the error is \(\sum_j \sigma_{\delta,j}^2 k_{x,j}(x,x) \bar h_j^{2p_j(x)}\). It goes to 0 as \(h \to 0\), and it is monotone in the componentwise order.
  - The orders are learned (Section 2.5).

### 2.4 Noise

$$ e_a \sim N(0, S_a), \quad S_a = D_a R_a D_a $$

- Runs are independent. Inside a run, \(R_a\) correlates the outputs.
  - [lit] Matrix-normal emulators count one run as one correlated vector: Overstall & Woods 1506.04489 eq 2. We keep \(R\) separate from the signal correlation.
- \(D_a\) is the diagonal of the noise sd: \(\log s^2(x, h) = m_s + \sum_j b_{s,j} \bar h_j + \zeta(x, h)\), with \(\zeta\) a latent GP.
  - [lit] hetGP 1611.05902 eq 14 smooths latent log-variances, so it works with few replicates.
- The trends \(b_{s,j}\) have symmetric priors, so the spread may grow or shrink as \(h \to 0\).
  - [lit] Analogous evidence: in an LES study the Lyapunov growth rate increases on finer meshes (1801.03046 §4.2).
- A replicate is a run with a small perturbation (for example a shifted start of the measurement). Whether that spread equals the physical spread is assumption A7.

### 2.5 Priors (all proper) and inference

| Parameter | Prior | Source |
|---|---|---|
| \(\beta\) (mean of \(\mu\)) | flat if no link and the basis matrix has full column rank; else Gaussian, centred on a physical estimate | improper posterior otherwise [inference] |
| \(c_{0j}\), \(c_{1j}\) | Gaussian with sd \(S_c\): the plausible size of the coarsest-level error relative to the converged value, in \(\Lambda\) units | proper, because \(c_{1j}\) multiplies \(\mu\) |
| \(\beta_s\) (module S8) | Gaussian | |
| \(p_{j0}\) | \(\log p_{j0} \sim N(0, 1)\) | wide; covers E&H's range |
| \(\pi_j\) (variation of \(p\) with \(x\)) | PC prior: \(P(\sigma_\pi > 0.3) = 0.05\), \(P(\ell < 0.1) = 0.05\) | shrinks to one shared order [assumption] |
| Matérn variance and range (\(r\), \(\delta\), \(\zeta\), \(R\)) | PC prior: \(P(\ell < 0.1) = 0.05\), \(P(\sigma > S) = 0.05\) | [lit] Fuglstad et al. 1503.00256 Thm 2.6, for an isotropic Matérn with \(d \le 3\); the paper says this limit "cannot be removed" (§2.3). Their range parameter equals \(2\ell\) in our Matérn form. Our use, one \(d = 1\) prior per input of an ARD kernel with up to 10 inputs, is a **heuristic, not a derived PC prior** [assumption]. G7 tests its influence. |
| \(\gamma\) (LB) | Uniform(0, 1) | |
| \(m_s\), \(b_{s,j}\) | Gaussian; \(b_{s,j}\) symmetric, with sd \(\log 4\) | [assumption] |

- The scales \(S\) are problem-specific, and they are the only problem-specific part (R7).
- A **known admissible range** (for example a positive rate) goes into the prior through a link, \(\mu = \psi(\tilde\mu)\). It does not go into a gate.
- **Sampling:** Gibbs.
  - Latent Gaussian fields use elliptical slice sampling. [lit] Murray, Adams & MacKay 1001.0175: it "has no free parameters" (§2.4). It needs a zero-mean Gaussian prior with a fixed covariance during the update (§2), so the mean levels (for example \(m_s\)) are sampled as separate variables (§2.5). It is efficient mainly when the prior dominates the likelihood (§2.5).
  - Their hyperparameters use the surrogate-data slice sampler. [lit] Murray & Adams 1006.0868 (abstract): it "requires little tuning while mixing well in both strong- and weak-data regimes". Its automatic form needs a likelihood that factorises over sites (§3.2). Given the other field, each of our two blocks does. The paper tested only fixed observation noise (§5), so its use with a latent log-variance field is [inference].
  - Censored values use data augmentation.
- **Convergence rule:** 4 chains. For every reported quantity, \(\hat R < 1.01\), where \(\hat R\) is the maximum of the rank-normalised split-\(\hat R\) and the folded split-\(\hat R\), and bulk-ESS > 400. [lit] Vehtari et al. 1903.08008 §2 and §4.1. We also require tail-ESS > 400 [assumption; the paper uses 400 for tail-ESS only in its examples].
- Why sample at all: [lit] ML-II "underestimate[s] prediction uncertainty" (1912.13440 §1).

### 2.6 Cost model

$$ \log_2 \kappa(u, h) = \kappa_0 + \sum_j \gamma_j \ell_j + \omega(u, h) + \eta $$

- \(\kappa\) is the cost per unit of run length, and \(c(a) = \kappa (T_{fixed} + T_{obs}) + c_0\).
- \(\omega\) is a GP. So the posterior can depart from the power law \(2^{\gamma \ell}\), which is only the prior mean. \(\gamma_j \sim N(3, 1)\) [assumption: 2D explicit time stepping, 4× cells and about 2× steps per level].
- [lit] Snoek 1206.2944 §3.2 puts a GP on log cost. Guinet 2011.11456 reports that simple cost models often predict better, so the prior mean carries most of the weight when data are few.
- **When the output is a time** (the run length needed is the unknown output): [lit] Hutter 1310.1947 Def. 1. The expected observed time is \(E[T_{obs}] = \int_0^T P_n(\tau > t) dt\).

### 2.7 Several candidate structures

- A structure is \(M = (k_h\) family, \(\Lambda\), set of levels used\()\).
- **Weights by stacking on an extrapolation score.**
  - Fit every \(M\) without the finest level, and score the held-out finest-level runs. This tests extrapolation in \(h\).
  - The scores are log densities **on the physical scale**, with the Jacobian of \(\Lambda\) included, so structures with different transforms are compared fairly.
  - [lit] Yao, Vehtari, Simpson & Gelman 1704.02030 eq 2.2 (log score). They recommend stacking because BMA "is flawed in the M-open setting".
  - **Our adaptation** [inference]: Yao et al. define and justify eq 2.2 with leave-one-out densities (§2.1–2.2). They do not discuss extrapolation hold-outs. We score a held-out finer level instead, because that task (predict one level finer than the data) is the closest available match to the real task (predict \(h = 0\)). This needs its own validation; assumption A13 records it.
  - Yao §2.3 warns that LOO "has large variance when the sample size is small". Our held-out sets are small, so the weights are noisy (see the minimum sample below).
  - A Dirichlet(2, ..., 2) penalty regularises the weights. [lit] Their §4.1 suggests "a strong prior … to the weights"; the specific Dirichlet(2, ..., 2) is our choice [assumption].
- **Cross-fitting:** the weights come from half of the held-out sites, and the calibration (Gate G1) from the other half; then swap.
- **Minimum sample:** each half needs at least 6 held-out sites [assumption]. With fewer, the weights stay equal, and G1 is reported as "not testable", never as "passed". Today the finest level (L10) has 1 run, so this rule applies.
- With only 2 levels left, \(p\) is not identifiable there. Then the weights stay equal and are flagged.

### 2.8 Prediction on the physical scale

- Pool the posterior draws over the structures with their weights, and back-transform them with \(\Lambda^{-1}\).
- \(m_y\) is the median. \(\sigma_y = (q_{84} - q_{16})/2\).
- Quantiles exist for every transform. For example, if a rate is Gaussian, the time \(1/\)rate has no finite moments.
- The Gaussian of the requirement is \(N(m_y, \sigma_y^2)\). Gate G6 checks its shape.
- \(\sigma_{epi}\) is \(\sigma_y\) of \(f(x)\).
- \(\sigma_{fid}(x,h) = (E[(f(x,h) - f(x))^2 \mid D])^{1/2}\), computed from the same draws.
- **Fidelity envelope**, in \(\Lambda\) units (relative error for \(\Lambda = \log\)):

$$ \sigma_{env}^2(x, h) = E_{\vartheta, c, \mu \mid D} [ \sum_j \bar h_j^{2 p_j(x)} ( (c_{0j} + c_{1j} \mu(x))^2 + \sigma_{\delta,j}^2 k_{x,j}(x, x) ) ] $$

  - The error of level \(h\) is \(\rho_0 + (\rho_1 - 1)\mu + \delta\). The envelope includes both its power-law part and its random part, so it does not understate the error when the mean explains most of it.
  - It is a sum of per-component squares, so there are no cross terms that could cancel. Each posterior draw is monotone in the componentwise order, so the average is monotone too. It is 0 at \(h = 0\), and it is finite because the posterior of \(c\) is proper.
- **Aleatoric spread on the physical scale:** \(s_0(x)\) is the quantile half-width \((q_{84} - q_{16})/2\) of \(\Lambda^{-1}(\Lambda(m_y) + e)\), with \(e \sim N(0, s^2(x, 0))\), over the posterior draws of \(s\). It is reported as "not identified" when replicates exist at fewer than 3 levels, or when the posterior sd of the trend \(b_s\) is more than half its prior sd.
- \(\sigma_{tot}\) is the quantile spread of draws of \(\Lambda^{-1}\)(target + noise) at \(h = 0\). This is the only definition; it is not \(\sigma_{epi}\) and \(s_0\) added in quadrature, because they are on scales that do not add.

---

## 3. Algorithm

**Step 0, set the problem.** Choose \(\Sigma\), \(\varepsilon\), \(C\) and the batch size \(q\). Choose the candidate transforms, an external cheap source if one exists (module S8), the prior scales, and the optional modules of Section 4.

**Step 1, initial design.**
- **Cheapest level:** a space-filling design, as large as is useful. This is Yi's abundant LF, and all of it enters the likelihood.
- **Prerequisite: the A12 test,** before the first real campaign. Step 5 values probes at unprobed finer levels correctly only if A12 holds.
  - Synthetic truths: 20 random draws of \(f(x, h) = f_0(x) + a(x) \bar h^{p}\), with \(x\) in 2D and \(p\) drawn from [0.7, 2.5]. There is data at \(\bar h = 1, 1/2, 1/4\), and candidates also at \(\bar h = 1/8\) and \(1/16\), with cost \(\propto 2^{3\ell}\).
  - Oracle: the value of each of the top 10 candidates, computed by full MCMC refits on its fantasy outcomes.
  - Pass rule: the acquisition's top candidate is in the oracle's top 3 in at least 16 of 20 cases [assumption].
  - If it fails: replace the reweighting by short MCMC refits for the 5 best candidates of each step. This costs more, and the cost must be measured.
- **Ladder:** nested Sobol designs in \(u\) at the levels above it. [lit] Stacking designs eq 9 gives starting sizes.
- **Replicates:** at 3 or more levels at 3 or more sites.
- This minimum is necessary, not sufficient [inference].

**Step 2, fit.** Run MCMC for every candidate \(M\) (Section 2.5), then compute the stacking weights (Section 2.7).

**Step 3, gates (diagnostics and calibration).**

| Gate | Test | Pass rule |
|---|---|---|
| G0 order | posterior of \(p_{j0}\); observed increment ratios | \(P(p_{j0} > 0.5) \ge 0.9\) **and** the posterior sd of \(\log p_{j0}\) is at most 0.5 (half its prior sd), so the prior cannot pass G0 alone (the prior already gives \(P(p > 0.5) = 0.76\)). When \(P(\sigma_\pi > 0.05 \mid D) \ge 0.9\) (the order varies with \(x\)), it also needs \(P(p_j(x) > 0.5) \ge 0.9\) at every \(x \in \Sigma_N\). **A failure does not stop the method.** It means the limit is not yet identified, so \(\sigma_{epi}\) is large, and Step 5 then prefers finer levels if A12 holds. [lit] E&H treat \(p < 0.5\) as anomalous. |
| G1 level hold-out | cross-fitted prediction of the finest level | 95% coverage within binomial limits; no sign bias; \(C_{LOO}\) near 1 [lit: Bachoc 1301.4320 eq 6; Oliver 1311.0828 eq 17] |
| G2 block LOO | leave out whole runs | z ~ N(0, 1); U statistic [lit: Overstall & Woods eq 10] |
| G3 noise | replicate spread against \(s^2\) | chi-square, p > 0.05 |
| G4 pre-asymptotic | refit without the coarsest level | the target moves less than \(\sigma_{epi}\); else remove that level |
| G5 structure | monotone coordinates stay monotone | no violation on \(\Sigma_N\) |
| G6 shape | 2.5/97.5% quantiles against \(m_y \pm 1.96\sigma_y\) | difference < 10% of \(\sigma_y\) |
| G7 prior | halve and double each prior scale (importance reweighting; refit if ESS < 400) | \(m_y\) moves < \(0.5\sigma_{epi}\) and \(\sigma_{epi}\) changes < 20%; else the output is flagged "prior-dominated" |

- **Repair:** at most 2 cycles per gate: remove a pre-asymptotic level (G4), recompute the weights, or buy a finer probe where \(|z|\) is largest.
- If a gate still fails, the output is labelled "uncalibrated". It is never silently reported as calibrated.
- The thresholds are [assumption].

**Step 4, stop or forecast.**
- **Success:** P1 holds and the gates pass.
- **Budget used:** report P2.
- **Forecast** (bounded cost):
  - Simulate the greedy policy of Step 5 with 20 fixed posterior samples (spread over the structures by their weights), for at most 50 steps or until success. Use variance-only updates and the expected reached set, with no fantasies and no reweighting.
  - Judge success on \(\sigma_{epi}\) of the **pooled** 20 samples, so the spread between samples and between structures (for example the 2–4× kernel spread of F6) stays in. Variance-only updates cannot shrink that spread, so the forecast is conservative about success.
  - Cost: about \(2 \times 10^8\) flops per sample and step (120 candidates, rank-10 updates at 400 data), so about \(2 \times 10^{11}\) flops for 20 samples and 50 steps.
  - **Known limit:** variance-only updates cannot shrink the spread between structures. When that spread exceeds \(\varepsilon\), the forecast says "P1 infeasible" from the first step and carries no information. A better forecast would also simulate how finer-level fantasies change the stacking weights. That is future work.
  - If the probability of exceeding \(C\) is above 0.5, tell the user the forecast. This is a report, not a stop: the method continues with the P2 criterion of Step 5.

**Step 5, acquisition: maximum uncertainty reduction per cost.**

$$ a^\star = \arg\max_{a \in A} \; (H_n - E J_n(a)) / E c(a), \qquad H_n = \sum_{x \in \Sigma_N} (\sigma_{epi}^2(x)/\varepsilon^2(x) - 1)_+ $$

- [lit] The ratio is MR-SUR (Stroh et al. 2007.13553 eq 11). The hinge form of \(H_n\) is our change [inference]. It is zero exactly when P1 holds.
- **P2 criterion.** When the forecast says P1 is infeasible, \(H_n\) becomes a soft maximum, \(\beta^{-1} \log \sum_{x \in \Sigma_N} \exp(\beta \sigma_{epi}^2/\varepsilon^2)\) with \(\beta = 20\) [inference]. This targets the max of P2 rather than a sum.
- **Budget enforcement.**
  - Every candidate gets a cap \(\bar c(a)\), the 0.95 quantile of its cost posterior. A candidate is admissible only if \(\bar c(a) \le C_{rem}\), where \(C_{rem}\) is the budget minus the spent cost and minus the caps of all pending runs (each pending run is reserved at its cap).
  - A run that reaches its cap is stopped. Its cost is counted, and it is a right-censored cost datum: the cost model (Section 2.6) uses a censored (Tobit) likelihood for it. The outputs it reached are kept; outputs not reached are censored (module S1). In the generic case with one output, a capped run gives only the cost datum. With module S2, it can be extended later if budget remains.
  - So \(\sum c \le C\) holds by construction.
- **Not optimal.** The policy is greedy, with a one-step look-ahead. It does not claim to minimise \(\sum c\) (P1) or to reach the P2 optimum. It is a heuristic in the class of MR-SUR.
- **Candidates:** every \(u\) in a candidate set, **at every resolution, including resolutions finer than any run so far**, and every optional action of Section 4.
  - The cost model prices each candidate, and the acquisition chooses. So refinement happens when it is the best use of the budget (Fig. E, curve 4).
- \(J_n(a)\) is a Monte Carlo average over 16–256 fantasies. Each fantasy:
  - draws the hyperparameters and the latent noise;
  - draws the outputs, including which nested outputs are reached and which are censored;
  - conditions each posterior draw's GP exactly on the fantasy outputs (Gaussian conditioning);
  - reweights the hyperparameter draws by the fantasy likelihood. If the effective sample size of the weights falls below 50 [assumption], it keeps the current weights for that candidate. That is conservative: it ignores the value of reducing the uncertainty of \(p\). This happens most for probes at unprobed finer levels;
  - keeps the stacking weights fixed [assumption], recomputes \(H\), and clips negative gains (Monte Carlo noise in a quantile spread) at 0.
  - [lit] Snoek 1206.2944 §3.3 uses Monte Carlo over pending outcomes.
  - The number of fantasies is doubled until the Monte Carlo error is below 10% of the gap between the two best candidates.
- **Batches:** choose greedily, and condition on the pending runs. [lit] Takeno 1901.08275 eq 7–8; Kathuria 1611.04088.

**Step 6, run the batch, add the data, and go to Step 2.**

**Computational cost** [inference, order of magnitude]: for \(n\) up to about 2,000 data in the likelihood, about 120 candidates, 16 fantasies and 200 posterior samples, rank-10 updates cost about \(4 \times 10^7\) flops each, so one acquisition costs about \(10^{13}\) flops: hours on one node. This is small against one L9 run. MCMC time must be measured.

**Output:** for each \(x \in \Sigma\): \(m_y\), \(\sigma_{epi}\), \(s_0\) (or "not identified"), and \(\sigma_{tot}\); \(\sigma_{env}(x, h)\) and \(\sigma_{fid}(x, h)\) for the levels used; the posteriors of \(p_j\) and \(c\); the stacking weights; the gate table; the spent cost and the allocation per level.

---

## 4. Optional modules: assumptions that use the structure of a problem

Each module is generic. The core (Sections 1–3) works without any of them. A module states an assumption, how it changes the algorithm, and how to check it.

| # | Structure | Assumption | Change to the algorithm | Check | Bioreactor |
|---|---|---|---|---|---|
| S1 | Nested outputs | One run returns the QoI at all values of an output coordinate up to a stopping value (a time-to-threshold curve, a load-displacement curve). | \(O(a)\) depends on the outputs; \(R\) correlates them; outputs not reached by \(T\) are censored. [lit] taKG 1903.04703 §2.3: shorter outputs come free. | G2 on whole runs | \(\Delta t_{mix}(\chi)\) for all χ ≤ the χ reached |
| S2 | Pause and resume | A run can continue from a checkpoint. | Action extend(\(j\), \(T \to T'\)) at incremental cost. [lit] Freeze-thaw 1406.3896; taKG §2.6. | restart parity test | chain/checkpoint tools exist |
| S3 | Warm start across levels | A finer run can start from a coarser state and lose less spin-up time. | Lower \(T_{fixed}\) for that action. | the result must not depend on the start (test A11) | the L10 runs were warm-started; test pending |
| S4 | Several resolution components | Parts of the model can be resolved separately. | \(h\) is a vector (Section 2.3); each part has its own order \(p_j\). | G0 per component | flow grid and scalar grid: a passive tracer on a finer grid, advected by a stored, converged, periodic flow |
| S5 | Several QoIs per run | One run gives many QoIs. | The criterion sums over QoIs, with weights \(1/\varepsilon_q^2\). [lit] Giles 1304.5472 §7.4. | — | kLa, \(\Delta t_{mix}\) and shear from one run |
| S6 | Known admissible range | The converged value must be in a known range. | a link in the prior (Section 2.5) | G7 | rates > 0 |
| S7 | Monotone coordinate | \(f\) is monotone in one coordinate. | Gate G5; enforce only if it is violated. [lit] López-Lopera 1901.04827 eq 8. | G5 | \(\Delta t_{mix}\) increases with χ |
| S8 | Cheap external source | A different, cheaper model (not a resolution of this code) correlates with the QoI. | \(\mu = \beta_s f_s + r\), with \(f_s\) a GP and its data in the likelihood (Section 2.2). | the posterior of \(\beta_s\) is away from 0 | none yet; L6 is a resolution level, so it is in the ladder |

---

## 5. Core assumptions

| # | Assumption | Status | Test |
|---|---|---|---|
| A1 | \(\Lambda(y)\) is Gaussian | untested | G6; Q–Q plots of G2 |
| A1b | The QoI types the method is meant for (mean or max shear, kLa, mixing time) follow A3 | tested here for \(\Delta t\) (fails on L6–L9) and kLa (plausible) only; untested for shear | G0 per QoI |
| A2 | The code converges, \(\delta(x, 0) = 0\) | fundamental; the data can show convergence only within the levels run | G0 |
| A3 | Leading power law in \(h\), with an order that is shared over \(X\) unless the data need \(p(x)\) | [assumption]. Bect Prop. 3 supports the power law for one QoI at one \(x\); nothing supports sharing across \(x\). Fig. D shows R from below 0 to above 10 across rpm. | G0; posterior of \(\sigma_\pi\) (Section 2.3) |
| A4 | The \(k_h\) family | stacking weights | G1 |
| A5 | Separable \(k_x k_h\) | untested | G2 residuals against \(x\) |
| A6 | Noise trend in \(h\) of either sign | [lit] analogy only | replicates |
| A7 | A perturbed replicate represents the physical spread | untested | experimental replicates |
| A8 | Stationary kernels in \(x\) | untested | G2 |
| A9 | (with S8) the external source is informative | posterior of \(\beta_s\) | reported |
| A10 | The coarsest ladder level is in the asymptotic range | per problem | G4 |
| A11 | (with S3) The warm start does not change the QoI | per problem | a direct comparison |
| A12 | Fantasy reweighting values the probes that reduce the uncertainty of \(p\) | [inference] | synthetic test |
| A13 | Stacking weights for one-level extrapolation transfer to \(h = 0\) | untestable without a finer level | a finer level |

---

## 6. Worked example: mixing time in a 2D rocking bioreactor

**Setting.**
- \(\Delta t_{mix}((\chi, \theta_{max}, rpm), h)\), with \(\chi \in [0.5, 0.95]\), \(\theta_{max} \in [2, 9]\)°, and rpm \(\in [15, 37.5]\).
- \(u = (\theta_{max}, rpm)\) and \(v = \chi\) (module S1).
- \(\varepsilon_{rel} = 0.10\).
- Budget: one week of the `priority` QOS, which is 312 cores, so \(C \le 52{,}000\) core-h.
- The physical system is taken to be 2D.

**Evidence so far:**

| # | Finding | Evidence |
|---|---|---|
| F1 | L6–L9 do not show convergence of \(\Delta t\). The increment ratio R must be about \(2^p > 1\) for convergence. Measured R: 0.2–0.4 for \(\Delta t\); 0.6–1.2 for \(\log \Delta t\); 1.5–3.4 for \(1/\Delta t\), but the Richardson limit of the rate is negative at 9 of 22 triplets. The limit is therefore not yet identified. | [data] Fig. D, `scripts/diag_observed_order.py` |
| F2 | kLa (no transform) converges plausibly: monotone at 90–100% of rpm on L6–L8 and 56–60% on L7–L9. | [data] `observed_order_kla.png` |
| F3 | No fixed kink in χ(t): the late/early decay-rate ratio is 0.36–3.65, and its direction changes with rpm and with level. | [data] Fig. A |
| F4 | Cost per level: ×7, ×9, ×21 per simulated second. Per probe to \(\Delta t_{0.95}\): ~0.015, 0.25, 6 and 300 core-h at L6–L9. The 80-cycle spin-up is 40–55% of an L9 run. | [data] Fig. C |
| F5 | L6 loses 7–10% of the tracer mass (L9: under 0.4%). | [data] diary 2026-10-03 |
| F6 | Held-out L9 from L6–L8: in log space all z-scores are positive (systematic). The kernel choice changes \(f(0)\) by 2–4×. | [data] Fig. B |
| F7 | Kim et al. ran uniform L10 (\(n_L = 2^{10}\), 0.24 mm, Main.tex:432), with L11 as the convergence reference in their appendix (user, 2026-10-03). The `MAXLEVEL = 9` in the imported driver is a default that our runs do not use (the level comes from `params.json`). Still open: Kim reports 120 core-h for one case, about 30× below our measured L10 cost. | `experiments/kimetal2024/Main.tex:432`, `:717` |
| F8 | **A short spin-up changes the QoIs** (test of A11; L8, 32.5 rpm). With release after 33 cycles: kLa +26% and \(\Delta t_{0.95}\) +12%. With release after 50 cycles: \(\Delta t_{0.75}\) +23%. These are against releases after 80–85 cycles, whose spread is under 1%. A warm start therefore must not shorten the settling. The only L10 point was released at absolute cycle 80, but through two checkpoint restarts. The restart at cycle 47 (`fig9_l10_seg1`) had no `*_prev` settings, which inject a seam transient (diary 2026-10-01). A test now repeats that restart history at L8, with and without `*_prev` (`scripts/submit_restart_test_l8.py`). Until it reports, the L10 point is not used. | [data] Fig. F, `scripts/plot_spinup_test.py` |
| F9 | **The run-to-run spread depends strongly on \(x\), and on \(h\) in no fixed direction** (releases 80/82/85, n = 3 per point). CV of \(\Delta t_{0.95}\), L6/L7/L8: 17/1.8/2.0% at 17.5 rpm; 1.2/0.15/1.9% at 25 rpm; 0.25/0.45/0.67% at 32.5 rpm. So a single flat noise value is wrong somewhere, the heteroscedastic noise model of Section 2.4 is needed, and the symmetric prior on the trend in \(h\) (A6) is supported. | [data] `scripts/diag_replicates.py`, `replicate_cv.png` |

![Fig. D](../experiments/multifidelity/observed_order_dtmix.png)

*Fig. D. Increment ratio R per rpm, for L6–L8 (circles) and L7–L9 (squares). Rows: \(\Delta t\), \(\log \Delta t\), \(1/\Delta t\). Dashed: R = 1, no convergence. Dotted: R = 2.8, order 1.5.*

![Fig. B](../experiments/multifidelity/h_kernel_test_chi95.png)

*Fig. B. \(\Delta t_{0.95}\) against cell size, per rpm: TWY2 and Brownian posteriors on L6–L9, with 95% bands. The cross is the single L10 run.*

![Fig. C](../experiments/multifidelity/cost_per_level.png)

*Fig. C. Measured cost. Left: core-s per simulated second. Right: core-h to observe \(\Delta t_{0.95}\) after release.*

![Fig. A](../experiments/multifidelity/chi_decay_rate_shape.png)

*Fig. A. The period-averaged decay rate \(\lambda(\chi)/\lambda(0.5)\) for L6–L9.*

**What the method will do with this** [inference]:
- G0 fails today, so \(\sigma_{epi}\) at \(h = 0\) is large.
- If assumption A12 holds (its synthetic test is a prerequisite), the acquisition will then value L10 probes, which are the only way to identify the limit.
- The Step 4 forecast will be uninformative here: the kernel spread (2–4×, F6) is far above \(\varepsilon_{rel} = 0.10\), so it will report "P1 infeasible" at once and the P2 criterion will start. That is the known limit stated in Step 4, not a finding about the problem.
- An L10 probe to \(\Delta t_{0.95}\) costs about 3,600 core-h at 32.5 rpm (spin-up included) and much more at low rpm. So within 52,000 core-h the method can buy roughly 5–14 of them, mostly at high rpm.
- [inference] With a Péclet number of about \(10^7\) (\(D = 0.44 \times 10^{-9}\) m²/s, `src/BioReactor.c:201`), the Batchelor scale is about 50 µm, against 0.49 mm at L9. If the asymptotic range starts only near that scale (L12–L13), the method will end in P2 with an honest, large \(\sigma_{epi}\). If L10 already shows convergence, it can succeed.
- Module S4 (a scalar grid finer than the flow grid on a stored converged periodic flow) could make the fine scalar levels much cheaper. It fits the general framework as a second resolution component. It needs a replay solver, which is in BACKLOG and is untested.

![Fig. F](../experiments/multifidelity/spinup_test_l8_32p5.png)

*Fig. F. Change of each QoI against the mean of the releases after 80, 82 and 85 cycles, for L8 at 32.5 rpm.*

**Runs** (25 jobs, submitted 2026-10-03, same binaries and protocol as the ladder):
- **Replicates:** tracer release at cycles 82 and 85 (and 80 at L7), at L6–L8 × 17.5/25/32.5 rpm and L9 × 25 rpm. They give \(s(x,h)\), its trend in \(h\), and \(R\). L6 and most of L7/L8 are done (F9); L9 is running.
- **Spin-up test** (A11): done, and it refutes A11 (F8).

---

## 7. Open questions

- **Q-a. Physical prior values, one set per QoI.** Today some come from looking at L6–L9, which uses the data twice; G7 tests their influence. The values below are what I need. A rough range from your physical judgement or the literature is enough; a guess with a stated source is better than my data-driven value.

| # | Quantity | Question | My current value (source) |
|---|---|---|---|
| 1 | \(S_\mu\) | Over \(\Sigma\), by what factor can the converged mixing time vary between its smallest and largest values (95% sure)? | a factor \(e^{2}\approx 7\) [assumption] |
| 2 | \(S_c\) | By what factor can the coarsest level used (L6) be wrong against the converged value (95% sure)? | a factor 50 (set after seeing L6–L9) |
| 3 | \(p\) | What convergence order do you expect for this QoI from the schemes (advection, VOF, embedded boundary)? | \(\log p \sim N(0, 1)\): median 1, 95% in [0.14, 7] [assumption] |
| 4 | \(b_{phys}\) | Only for the rate link: a physical estimate of the mixing time in rocking cycles, from experiments | 30 cycles at 25 rpm (no source) |
| 5 | noise | Run-to-run spread | **no longer needed**: measured (F9) |
- **Q-b.** Should module S4 (two resolution components with a replayed flow) be developed now, or after L10 probes show whether the scalar converges?
- **Q-c, answered 2026-10-03:** fully Bayesian; no deterministic plug-in (Section 2.2).
- **Q-d, answered 2026-10-03:** the median target is accepted.
