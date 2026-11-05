---
theme: default
title: Multi-fidelity Bayesian optimization of mixing in two-phase rocking bioreactors
aspectRatio: 16/9
canvasWidth: 1024
mdc: true
slideNumber: true
colorSchema: light
transition: fade
---

<video src="/title_v2_light.mp4" autoplay loop muted playsinline class="absolute inset-0 w-full h-full object-cover" />

<div class="absolute left-14 top-12 right-14">
  <div style="font-size: 40px; line-height: 1.15; font-weight: 300; letter-spacing: -0.5px; color: #14202e">Multi-fidelity Bayesian optimization of mixing in two-phase rocking bioreactors</div>
  <div class="mt-5" style="font-size: 22px; font-weight: 200; color: #14202e"><b style="font-weight: 300">Elvis Aguero</b> · Minki Kim · Radu Cimpeanu · Miguel Bessa · Daniel Harris</div>
  <div class="mt-1" style="font-size: 17px; font-weight: 200; color: #7a8797">Brown University · Kookmin University · University of Warwick</div>
</div>

---
transition: fade
clicks: 1
---

<div class="slide-title">cultivated meat needs large, gentle bioreactors</div>
<div class="mot-all">
<div class="mot-figs">
<div class="mot-pcol">
<img src="/cm_process.png" class="mot-proc" />
<div class="cap">(Foodingredients 2021; Post et al. 2020; Guan et al. 2021)</div>
</div>
<div class="mot-vcol">
<video src="/celltainer_rock.mp4" autoplay loop muted playsinline class="mot-vid"></video>
<div class="cap">CELL-tainer Biotech BV, youtube.com/watch?v=jC-eq9N16o4</div>
</div>
</div>

<div class="mot-band">
<div class="mot-cost">
<span class="big">$25–35</span><span class="u">/kg</span>
<span class="arrow">→</span>
<span class="big">$9</span><span class="u">/kg</span>
<span class="cap">Negulescu 2023</span>
</div>
<div class="mot-cfd" v-click="1">
<span class="name">one resolved condition costs</span>
<span class="big2">≈ 4,000</span><span class="u">CPU-h</span>
</div>
</div>

<div class="mot-eq">

$$
\underbrace{\tau < \tau_{\text{damage}}}_{\text{1–100 Pa}}
\qquad
\underbrace{k_La\,(C_{\text{sat}} - C) \;\ge\; \mathrm{OUR}}_{\text{oxygen}}
\qquad
\underbrace{t_{\text{mix}} < 1/k_La}_{\text{mixing}}
$$

</div>

<div class="key">

$\tau$ shear stress · $\tau_{d}$ cell-damage threshold · $k_La$ oxygen transfer coefficient · $C$ dissolved oxygen · $C_{\text{sat}}$ saturation value · OUR oxygen uptake rate of the cells · $t_{\text{mix}}$ mixing time

</div>
</div>

<style>
.mot-all { position: absolute; left: 32px; right: 32px; top: 56px; display: flex; flex-direction: column; }
.mot-figs { display: flex; gap: 24px; align-items: flex-start; }
.mot-pcol { width: 476px; flex: none; }
.mot-vcol { width: 460px; flex: none; }
.mot-proc { display: block; width: 476px; margin: 0 !important; }
.mot-vid { display: block; width: 460px; margin: 0; }
.mot-figs .cap { margin-top: 4px; }
.mot-band { display: flex; justify-content: space-between; align-items: center; margin-top: 8px; }
.mot-cost { display: flex; align-items: baseline; gap: 12px; white-space: nowrap; }
.mot-cost .big { font-size: 36px; font-weight: 300; color: var(--ink); }
.mot-cost .u { font-size: 20px; color: var(--mute); }
.mot-cost .arrow { font-size: 34px; color: var(--mute); }
.mot-cost .cap { align-self: center; margin-left: 4px; white-space: nowrap; }
.mot-cfd { display: flex; align-items: baseline; gap: 10px; white-space: nowrap; }
.mot-cfd .big2 { font-size: 36px; font-weight: 300; color: var(--ink); }
.mot-cfd .u { font-size: 20px; color: var(--mute); }
.mot-eq { margin-top: 13px; text-align: center; }
.mot-eq .katex-display { margin: 0; }
.mot-eq .katex { font-size: 1.0em; }
.key { margin-top: -12px; } /* trims the KaTeX display strut: visible gap to the equations is about 24 px */
</style>

<!--
Cultivated meat is still far from the price of beef. Negulescu and co-workers estimate 25 to 35 dollars per kilo in stirred tanks of 42 to 211 cubic metres. To compete with beef, production must drop below 9 dollars per kilo. Media dominate the cost, but the bioreactor sets the cell density and the scale.

The bioreactor must do three things at once. The shear stress must stay below the damage range, 1 to 100 pascal for mammalian cells. The oxygen transfer must supply the oxygen uptake of the cells. The mixing time must be shorter than one over kLa, so the oxygen spreads before the cells use it. These pull against each other: more agitation helps oxygen and mixing, but raises shear. Humbird 2021: "shear-induced cell damage will all limit practical bioreactor volume and attainable cell density".

[click] Simulations can map this trade-off before experiments. But one resolved condition costs about 4,000 CPU hours at our finest grid. Zhan 2019 say this directly for their wave bag: "highly time consuming and not practically feasible for a large number of operating conditions".

The video shows a rocking bioreactor: a single-use bag on a platform that rocks by a few degrees. The free surface makes waves. The waves mix the medium and bring oxygen in, with no impeller and no bubbles. Mosa Meat uses this CELL-tainer system. Note: the CELL-tainer also translates; our model rocks only.
-->

---
transition: fade
---

<div class="slide-title">earlier work</div>
<div style="position:absolute; left:0; top:44px; width:100%; height:100%">
<img src="/s3_kim.png" class="s3-img" style="left:16px; top:58px; width:269px; height:85px; object-fit:cover" />
<div class="s3-lab" style="left:16px; top:200px; width:269px"><div class="au">Kim et al. 2025</div><div class="ct">resolved Basilisk model</div></div>

<img src="/s3_seidel.png" class="s3-img" style="left:297px; top:10px; width:310px; height:175px" />
<div class="s3-lab" style="left:297px; top:200px; width:319px"><div class="au">Seidel et al. 2022</div><div class="ct">kLa and mixing maps</div></div>

<img src="/s3_zhan.png" class="s3-img" style="left:628px; top:10px; width:380px; height:180px" />
<div class="s3-lab" style="left:628px; top:200px; width:380px"><div class="au">Zhan et al. 2019</div><div class="ct">first CFD of a wave bag</div></div>

<img src="/s3_cantarero.png" class="s3-img" style="left:40px; top:278px; width:255px; height:150px" />
<div class="s3-lab" style="left:40px; top:438px; width:360px"><div class="au">Cantarero-Rivera et al. 2024</div><div class="ct">coarse-to-fine correction</div></div>

<img src="/s3_martens.png" class="s3-img" style="left:440px; top:280px; width:530px; height:153px" />
<div class="s3-lab" style="left:440px; top:436px; width:530px"><div class="au">Martens et al. 2025</div><div class="ct">multi-fidelity GP, bioprocess</div></div>

</div>
<div class="s3-close">Our goal: predict the flow over all operating conditions from cheap simulations, with error bars.</div>

<style>
.s3-img { position: absolute; display: block; }
.s3-lab { position: absolute; font-size: 16px; line-height: 1.25; color: var(--ink); }
.s3-lab .au { font-size: 16px; }
.s3-lab .ct { color: var(--mute); font-size: 16px; }
.s3-close { position: absolute; left: 16px; right: 16px; top: 530px; font-size: 18px; color: var(--ink); border-top: 2px solid #e8590c; padding-top: 6px; }
</style>

<!--
Five earlier works. Each one gave us one piece of the problem.

Zhan 2019 was the first resolved CFD of a 10 litre wave bag. They used a volume of fluid model in Fluent. They ran nine rocking speed and angle settings in 2D, and one 3D check (pages 58, 59, 61).

Seidel 2022 modelled the CELL-tainer, a bag with two degrees of motion. They tracked the motion with a camera and scanned the bag in 3D. They fitted response surfaces for kLa and mixing time (pages 2 to 4, 12, 13).

Kim 2025 is our solver base. It is a two-phase Basilisk solver in the frame of the bag, so the walls stay fixed. It transports oxygen across the interface. The kLa comes from the full oxygen field (pages 3 to 6, 9).

Cantarero-Rivera 2024 corrects a coarse CFD mesh of a spinner flask with a neural network. The model uses a steady RANS k-epsilon spinner flask, with two meshes of about 44,500 and 1.03 million elements. It is trained at one speed, 90 rpm. It lowers the velocity error by about 20 percent (abstract, page 140). CFD takes 2.5 hours; the coarse mesh and the network take 15 minutes (page 141). At 30 rpm the correction makes the error worse, by 99 percent in RMSE (Table 5). There are no error bars.

Martens 2025 uses a multi-fidelity Gaussian process for batch Bayesian optimization. The fidelities are three reactor scales of one CHO model. The study is simulation only. The cost cut of 93 percent is from one run against industrial design of experiments (section 3.3).

Our goal: predict the flow over all operating conditions from cheap simulations, with error bars.
-->

---
clicks: 4
---

<div class="slide-title">two-phase Navier–Stokes, in the frame of the bag</div>
<div class="mdl-left">

$$
\begin{aligned}
\underbrace{\rho\left(\partial_t\mathbf{u} + \mathbf{u}\cdot\nabla\mathbf{u}\right)}_{\text{inertia}}
&= \underbrace{-\nabla p}_{\text{pressure}}
+ \underbrace{\nabla\cdot(2\mu\mathbf{D})}_{\text{viscosity}}
+ \underbrace{\rho\,\mathbf{g}(t)}_{\text{gravity}}
+ \underbrace{\gamma\kappa\,\mathbf{n}\,\delta_s}_{\text{surface tension}} \\
&\quad \underbrace{{\color{#D55E00}-\,\rho\left(\dot{\boldsymbol\Omega}\times\mathbf{x} + 2\boldsymbol\Omega\times\mathbf{u} + \boldsymbol\Omega\times\boldsymbol\Omega\times\mathbf{x}\right)}}_{\color{#D55E00}\text{frame: angular acceleration, Coriolis, centrifugal}}
\end{aligned}
$$

<div class="key">

$\mathbf{u}$ velocity · $p$ pressure · $\rho$ density · $\mu$ viscosity · $\mathbf{D}$ rate of strain · $\mathbf{g}(t)$ gravity in the bag frame · $\gamma$ surface tension · $\kappa$ curvature · $\mathbf{n}$ interface normal · $\delta_s$ interface delta · $\boldsymbol\Omega$ rocking rate · $\mathbf{x}$ position

</div>

<div class="mdl-row">
<div>

$$\partial_t\alpha + \mathbf{u}\cdot\nabla\alpha = 0$$

<div class="cap">volume of fluid</div>
</div>
<div>

$$c_g = H\,c_\ell$$

<div class="cap">Henry's law at the interface</div>
</div>
</div>

<div class="key mdl-key2">

$\alpha$ liquid fraction · $c_g$, $c_\ell$ gas and liquid oxygen · $H$ Henry constant

</div>

<div class="mdl-grid">
<img src="/basilisk.png" style="height:40px; margin-bottom:0" />
<div>

$${\color{#0072B2}h_8 = \tfrac{L}{2^{8}}} \quad {\color{#009E73}h_9 = \tfrac{L}{2^{9}}} \quad {\color{#CC79A7}h_{10} = \tfrac{L}{2^{10}}}$$

<div class="cap">Basilisk, uniform grid, no refinement</div>

<div class="key">

$h_\ell$ cell size at level ℓ · $L$ tank length

</div>

</div>
</div>

</div>

<div class="mdl-right" v-click="1">
<div class="col-label">multi-fidelity</div>

$${\color{#CC79A7}f_H(x)} = \rho\,{\color{#0072B2}f_L(x)} + b + \delta(x), \qquad \delta \sim \mathcal{GP}(0,\, s^2 k_\ell) $$

<div class="mdl-fig">
  <img v-if="$clicks <= 1" src="/mf_stage1.png" />
  <img v-else-if="$clicks == 2" src="/mf_stage2.png" />
  <img v-else-if="$clicks == 3" src="/mf_stage3.png" />
  <img v-else src="/mf_stage4.png" />
</div>

<div class="key mdl-key3">

$f_H$ high-fidelity output · $f_L$ low-fidelity output · $\rho$ scale factor · $b$ offset · $\delta$ discrepancy · $s$, $k_\ell$ its amplitude and kernel · $\theta$ hyperparameters · $\sigma_n$ noise · $y_H$ high-fidelity data

</div>

<div v-click="3">

$$\theta^\star = \arg\min_\theta\, \left[-\log p(\mathbf{y}_H \mid \theta)\right], \qquad \theta = (\rho,\, b,\, \ell,\, s,\, \sigma_n)$$

<div class="cap">Yi, Cheng & Bessa 2026, Comput. Mech.</div>

</div>
</div>

<style>
.mdl-left { position: absolute; left: 16px; top: 62px; width: 540px; }
.mdl-right { position: absolute; right: 10px; top: 62px; width: 450px; }
.mdl-left .katex-display, .mdl-right .katex-display { margin: 6px 0; }
.mdl-left .katex-display .katex { font-size: 1.05em; }
.mdl-left .katex-display { margin: 14px 0 !important; }
.mdl-row { display: flex; justify-content: space-around; margin-top: 0; text-align: center; }
.mdl-grid { display: flex; flex-direction: row; align-items: center; gap: 14px; margin-top: 0; }
.mdl-grid > div { text-align: center; }
.mdl-grid img { height: 130px; }
.mdl-cost { margin-top: 4px; font-size: 17px; white-space: nowrap; color: var(--ink); }
.mdl-fig img, .mdl-fig video { width: 320px; display: block; }
.mdl-left .mdl-grid .katex-display { margin: 0 !important; }
.mdl-right .katex-display .katex { font-size: 0.86em; }
.mdl-right .katex-display { overflow: visible; }
.mdl-left .key { margin-top: 2px; }
.mdl-key2 { margin: 2px 0 0; }
.mdl-right .mdl-key3 { margin: 0 0 4px; }
</style>

<!--
The solver is Basilisk. We solve the two-phase Navier–Stokes equations in the frame of the bag, so the walls do not move. The price is the frame terms, in orange: the angular acceleration, Coriolis and centrifugal terms. The volume fraction alpha marks the liquid. At the interface the oxygen follows Henry's law. The mesh is a uniform grid with no refinement. L8, L9 and L10 mean a cell size of L over 2 to the 8, 9 and 10.

[click] Multi-fidelity. Blue dots: many cheap L8 runs. Pink dots: few expensive L10 runs. Measured cost. One L10 run costs about 65 core-hours per rocking cycle (median of the 13 training runs); one L8 run costs about 0.22. The cost slide uses about 670 core-hours per warm L10 run. So we can afford few L10 runs.

[click] A kernel ridge regression fits the cheap data: f_L.

[click] The expensive data are rho times f_L, plus b, plus a Gaussian process delta. We fit rho, b and the GP parameters by minimising the negative log-likelihood. Each pink curve is one step of the optimiser. The newest step is the strong line. On the right, the bold parameters are the ones that change at that step. At the end we get the prediction and its 95 % band. The data here are synthetic. Our data are in the results. We chain this: L8 to L9, then L9 to L10.
-->

---
clicks: 3
---

<div class="slide-title">what we measure</div>
<div class="g4">
<div class="cell">
<div class="head"><span class="name">shear stress</span><div class="eq4">

$\tau = \mu\left(\dfrac{\partial u}{\partial y} + \dfrac{\partial v}{\partial x}\right)$

</div></div>
<div class="vid">
  <SwapVideo loop="/tau_loop.mp4" rel="/oxy_rel.mp4" :on="$clicks >= 1" />
  <div class="cbar"><span class="bar bar-pink" v-if="$clicks < 1"></span><span class="bar bar-blue" v-else></span><span class="lab" v-if="$clicks < 1">shear stress <i>τ</i></span><span class="lab" v-else>oxygen <i>c</i></span></div>
</div>
</div>
<div class="cell" v-click="1">
<div class="head"><span class="name">oxygen transfer</span><div class="eq4">

$\dfrac{\mathrm{d}C^*}{\mathrm{d}t} = k_La\,(1 - C^*)$

</div></div>
<img src="/def4_kla.png" class="plot" />
<div class="note4">

$C^* = \dfrac{1}{V_\ell}\displaystyle\int_{V_\ell} \dfrac{c}{c_{\mathrm{sat}}}\,\mathrm{d}V$

</div>
</div>
<div class="cell" v-click="2">
<div class="head"><span class="name">dissipation rate</span><div class="eq4">

$\varepsilon = 2\mu\,\mathbf{S}\!:\!\mathbf{S}$

</div></div>
<div class="vid">
  <SwapVideo loop="/edr_loop.mp4" rel="/tracer_rel.mp4" :on="$clicks >= 3" />
  <div class="cbar"><span class="bar bar-orange" v-if="$clicks < 3"></span><span class="bar bar-blue" v-else></span><span class="lab" v-if="$clicks < 3">dissipation <i>ε</i></span><span class="lab" v-else>tracer <i>c</i></span></div>
</div>
</div>
<div class="cell" v-click="3">
<div class="head"><span class="name">mixing time</span><div class="eq4">

$\chi = 1 - \dfrac{\sigma^2}{\sigma_0^2}$

</div></div>
<img src="/def4_mix.png" class="plot" />
<div class="note4">

$\sigma^2 = \dfrac{1}{V_\ell}\displaystyle\int_{V_\ell} (c - \bar{c})^2\,\mathrm{d}V$

</div>
</div>
</div>

<style>
.g4 { position: absolute; inset: 70px 34px 20px 34px; display: grid; grid-template-columns: 1fr 1fr; grid-template-rows: 1fr 1fr; gap: 14px 40px; }
.cell { position: relative; display: flex; flex-direction: column; min-height: 0; }
.head { display: flex; align-items: center; gap: 18px; height: 62px; }
.head .name { font-size: 26px; }
.eq4 p, .note4 p { margin: 0; }
.eq4 .katex { font-size: 1.32em; color: var(--ink); }
.vid { position: relative; flex: 1; display: flex; flex-direction: column; justify-content: center; }
.plot { height: 178px; margin: 0 !important; align-self: flex-start; }
.note4 { position: absolute; left: 64px; top: 66px; }
.note4 .katex { font-size: 0.95em; color: var(--mute); }
.cbar { display: flex; align-items: center; gap: 8px; align-self: flex-end; margin-top: 8px; }
.cbar.shown { opacity: 1; }
.cbar .bar { width: 110px; height: 8px; border-radius: 4px; }
/* colours sampled from the video frames: light to dark, as in each video's colormap */
.cbar .bar-pink { background: linear-gradient(90deg, #f0d6e7, #eba9cb, #e05a91, #441e6d); }
.cbar .bar-orange { background: linear-gradient(90deg, #f7e9d7, #fac792, #f18843, #da6728, #68211e); }
/* after the click the videos show the release colormap (blue), sampled from oxy_rel and tracer_rel late frames */
.cbar .bar-blue { background: linear-gradient(90deg, #eef3f8, #b9dcee, #5fb7d6, #1f6fc0, #1d2a78); }
.cbar .lab { font-family: "Helvetica Neue", "Inter", sans-serif; font-weight: 200; font-size: 20px; color: var(--ink); }
</style>

<!--
Each quantity comes from one run: L8, 32.5 rpm, 7 degrees. The videos are in real time and in the lab frame.

Shear stress hurts the cells. The colour is the magnitude of the shear stress. We report the 95th percentile over the liquid.

[click] kLa is the speed of the oxygen transfer. At this click the oxygen starts to enter at the free surface (left video). C* is the mean dissolved oxygen in the liquid, divided by the saturation value. We fit the law over one rocking period at C* = 0.10. Here kLa = 47 per hour.

[click] The dissipation rate is the loss of kinetic energy in the liquid. S is the strain-rate tensor, one half of (grad u + grad u transpose). The colour is the log of epsilon.

[click] Mixing time. At this click the tracer is released in the top half of the liquid (left video). sigma squared is the variance of the tracer over the liquid; c bar is its mean. sigma_0 squared is the variance at the release. chi goes from 0 (not mixed) to 1 (fully mixed). We stop the clock at chi = 0.5: here 11.2 seconds.

Protocol: we settle the flow for 80 rocking cycles. Then we release the tracer in the top half of the liquid and start the oxygen at the same time.
-->


---
transition: fade
---

<img src="/s6_slide_r11.png" class="s6-bg" />

<div class="s6-head">Extending the multi-fidelity algorithm to three levels</div>

<div class="cap s6-cite">extends Yi, Cheng & Bessa 2026, Comput. Mech.</div>

<div class="s6-eq" style="top:56px">

$$f_8(x) = \mathrm{KRR}(\mathcal{D}_8)$$

</div>
<div class="s6-sub" style="top:128px">0.22 core-h/cycle</div>

<div class="s6-eq" style="top:165px">

$$f_9(x) = \beta_9 + \rho_9\, f_8(x) + \delta_9(x)$$

</div>
<div class="s6-sub" style="top:236px">5.76 core-h/cycle</div>

<div class="s6-eq" style="top:275px">

$$f_{10}(x) = \beta_{10} + \rho_{10}\, f_9(x) + \delta_{10}(x)$$

</div>
<div class="s6-sub" style="top:346px">64.9 core-h/cycle</div>

<div class="key s6-key">

$f_\ell$ prediction at level $\ell$ · KRR kernel ridge regression · $\mathcal{D}_8$ L8 training data · $\beta_\ell$ offset · $\rho_\ell$ scale factor · $\delta_\ell$ GP correction · GP Gaussian process · $w_g$ weight of grid point $g$ · $\lambda_g$ kernel length scale · $\eta_g$ noise level · $y$ training values

</div>
<div class="s6-strip">

$$\delta_\ell \sim \mathcal{GP}, \qquad w_g \propto p(\lambda_g)\, p(\eta_g)\, p(y \mid \lambda_g, \eta_g)$$

</div>

<style>
.s6-head { position: absolute; left: 32px; top: 12px; font-size: 26px; line-height: 1.5; color: var(--ink); }
.s6-bg { position: absolute; left: 0; top: 0; width: 1024px; height: 576px; }
.s6-eq { position: absolute; left: 32px; width: 450px; font-size: 30px; color: var(--ink); text-align: left; }
.s6-eq .katex { font-size: 1.0em; }
.s6-sub { position: absolute; left: 32px; width: 450px; font-size: 22px; color: var(--ink); }
.s6-key { position: absolute; left: 32px; width: 440px; top: 392px; }
.s6-strip { position: absolute; left: 0; right: 0; top: 508px; text-align: center; font-size: 26px; color: var(--ink); }
.s6-cite { position: absolute; left: 32px; top: 54px; }
</style>

<!--
Rows, top to bottom: L8, L9, L10. Each row has the equation for that level on the left, and the cropped tank on the right.
Row 1 (L8): f_8 is a kernel ridge fit on 40 cold L8 runs. Median cost 0.22 core-h per simulated cycle (IQR 0.16 to 0.26).
Row 2 (L9): the chain trains on 16 L9 runs. Six are measured warm runs (x8_l9) at 5.76 core-h per cycle (wall-clock logs; 11.8 simulated cycles each). The other 10 are cold fig9 theta-7 runs, costed as warm estimates: 5.76 core-h per cycle times 11.8 cycles, about 68 core-h per point, 682 core-h in total. The 16-run figure includes these 10.
Row 3 (L10): 13 warm runs. Median 64.9 core-h per cycle, steady rate excluding the first two cycles. Each run simulates about 6 cycles (x8 about 12). The logs do not separate start-up from steady running, so the rate is not corrected for it.
Bottom: the level below is the mean basis (beta, rho) of the level above, with a Gaussian process delta. The weights over the length scale lambda and noise eta are the prior times the marginal likelihood. Code: krr_lr_gpr_fullbayes.py lines 99 to 106 (basis), 126 (K_g), 137 (a_n, b_n), 140 (log marginal likelihood), 142 (weights), 169 to 171 (prediction). The level below enters as a plug-in mean (h3_band_check.py lines 23 to 27).
Tanks: cropped to the liquid plus a thin gas band, equal aspect. One shared colour scale: the 10th and 90th percentiles of |u| over the pooled liquid cells of the three levels, each level weighted equally (vmin 0.0228, vmax 0.2026, nondimensional speed). Saturated liquid cells above vmax: L8 8.4 %, L9 12.0 %, L10 7.2 %. Below vmin: L8 16.6 %, L9 5.1 %, L10 9.5 %. The upper limit moved from the 85th to the 90th percentile because the 85th left 19 % of L9 saturated.
-->

---
transition: fade
---

<div class="s8-head">accuracy against cost</div>
<img src="/s8_cost.png" class="s8-img" />

<style>
.s8-head { position: absolute; left: 32px; top: 14px; font-size: 26px; line-height: 1.5; padding-bottom: 6px; color: var(--ink); }
.s8-img { position: absolute; left: 42px; top: 46px; width: 940px; }
.s8-cap { position: absolute; left: 32px; right: 32px; top: 470px; font-size: 18px; line-height: 1.3; color: var(--ink); }
</style>

<!--
x axis: training cost in core-hours, log scale.
Plotted: the full-set chain only (solid line, filled markers), with all runs warm-started. Its cost is 40 L8 runs and 16 L9 runs (10 estimated warm, 6 measured warm) as the base, 3,170 core-h, then 667 core-h per L10 run. Each point is one number of L10 runs, n = 3, 4, 6, 8, 10, 12.
Also shown: H1 (2,043 core-h base, then 703 per L10 run) and L10 GP only (726 core-h per run).

Variant table (cheap chain subsets, not plotted). Errors are relative RMS, held-out. Errors at n = 3, 4, 6, 8, 10, 12 L10 runs.
| variant | n_L8 | n_L9 | base core-h | tau_95 error (%) | EDR error (%) |
|---|---|---|---|---|---|
| full set | 40 | 16 | 3,170 | 5.4, 5.8, 5.4, 5.2, 4.3, 2.3 | 8.7, 6.8, 8.9, 8.1, 6.6, 5.3 |
| warm9 | 40 | 6 | 2,488 | 9.0, 9.7, 6.4, 5.7, 5.1, 2.7 | 23.9, 43.5, 33.1, 23.9, 17.8, 11.9 |
| warm9, 20 L8 | 20 | 6 | 1,455 | 6.9, 6.9, 5.1, 5.1, 4.7, 3.4 | 48.9, 504, 58.9, 42.3, 28.4, 15.5 |
| warm9, 10 L8 | 10 | 6 | 918 | 5.6, 5.2, 4.9, 5.1, 4.0, 3.3 | 53.6, 1,080, 64.6, 43.5, 31.3, 18.9 |
Each cheap variant adds 703 core-h per L10 run.

Crossings, full-set chain only (core-h):
- tau_95 at 10 %: 5,170 core-h (n = 3, error 5.4 %).
- tau_95 at 20 %: 5,170 core-h (n = 3).
- tau_95 at 5 %: 9,837 core-h (n = 10, error 4.3 %). Not reached at n = 6 or 8.
- EDR at 10 %: 5,170 core-h (n = 3, error 8.7 %).
- EDR at 20 %: 5,170 core-h (n = 3).
- EDR at 5 %: not reached by n = 12 (5.3 % at n = 12).

Kim's code, measured on Oscar (job 8363941, N=1024, 64 ranks, steady t=0.1–0.4): ≈102 core-h per rocking cycle, against 26.3 for ours (v_l10_rpm32.5_th7); at N=256, 0.98 against 0.22.
-->

---
transition: fade
---

<div class="slide-title">2d campaign, baseline data</div>
<img src="/s9_heatmaps.png" class="s9-img" />

<style>
.s9-img { position: absolute; left: 62px; top: 46px; width: 900px; }
</style>

<!--
This slide shows the 2D campaign as a template. The data are the baseline runs only. No new run is included.

Rows are the three gated quantities: tau 95, the energy dissipation rate, and kLa at 10 percent oxygen. Columns are the predicted median on the left, and the relative uncertainty on the right.

The uncertainty is two sigma over the mean. The black line marks 0.30, the gate for a point to count as known.

The coloured dots are the runs. Blue is L8, green is L9 and pink is L10.

The energy dissipation rows use the time mean as a stand-in. The 95th percentile (edr_95) is not yet saved for the existing runs. Its backfill jobs are queued and not finished. The final version will use edr_95.
-->


---

<div class="slide-title">Conclusions</div>

---

<div class="slide-title">Questions?</div>

---
transition: fade
---

<div class="slide-title">backup</div>

---
transition: fade
---

<div class="slide-title">earlier work, by feature</div>
<img src="/lit_gap.png" class="lit-bk" />

<style>
.lit-bk { position: absolute; left: 60px; top: 50px; height: 470px; margin-top: 0 !important; }
</style>

<!--
What earlier work did. Zhan, Seidel and Kim simulate rocking bags with a resolved interface, each on one grid, with no error bars. Cantarero-Rivera corrects a coarse mesh with a neural network, in a spinner flask, again with no error bars, and it does not generalise across speeds. Kiviat computes all three quantities, but in stirred tanks with RANS, which they say under-predicts shear. Martens uses a multi-fidelity Gaussian process, but the fidelities are culture scales from an ODE model, not grid levels. We do all five: rocking bag, resolved two-phase flow, all three quantities, they correct one field with two meshes at one training speed, without error bars; we use three grid levels across the operating range, with error bars.
-->

---
transition: fade
---

<div class="s7-head">validation and held-out prediction</div>
<img src="/s7a_band.png" class="s7-a" />
<img src="/s7b_predict.png" class="s7-b" />

<style>
.s7-head { position: absolute; left: 32px; top: 14px; font-size: 26px; line-height: 1.5; color: var(--ink); }
.s7-a { position: absolute; left: 47px; top: 56px; width: 930px; }
.s7-b { position: absolute; left: 47px; top: 318px; width: 930px; }
</style>

<!--
Top: our L10 runs and Kim et al. 2025 at 7 degrees. The pink band is plus and minus two sigma of the aleatoric spread. Each point has its own sigma, from its own run: we resample the run's rocking cycles and recompute the three-cycle statistic. Sigma is 0.1 to 2.6 percent; it is largest at 22.5 to 25 rpm. So the band is thin. Our shear amplitude is 10 to 15 percent above Kim from 25 to 32.5 rpm, far outside the band: this gap is systematic, not noise. The dissipation agrees with Kim. Limits: each sigma comes from four cycles only. At 22.5 rpm (both statistics) and 32.5 rpm (dissipation) the cycles still drift, so part of that sigma is transient, not noise.

Bottom: the model predicts each run before it sees that run. The bars show the 95 % band.

Left: the 13 training runs, each left out once.

Right: six new runs at 5.5 degrees. Their predictions use no refit.

The EDR is the maximum over time, as Kim defines it. At the new angle, the time mean is more accurate than this maximum.
-->
