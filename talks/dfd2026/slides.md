---
theme: default
title: Multi-fidelity Bayesian optimization of mixing in two-phase rocking bioreactors
aspectRatio: 16/9
canvasWidth: 1024
mdc: true
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
clicks: 2
---

<div class="mot-left">
<div class="name">cultivated meat needs large, gentle bioreactors</div>

<div class="mot-cost">
<div><span class="big">$25–35</span><span class="u">/kg</span><div class="cap">stirred tanks, 42–211 m³<br/>Negulescu 2023</div></div>
<div class="arrow">→</div>
<div><span class="big">$9</span><span class="u">/kg</span><div class="cap">to compete with beef<br/>Negulescu 2023</div></div>
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

<div class="quote">"shear-induced cell damage will all limit practical bioreactor volume and attainable cell density"<span class="who"> — Humbird 2021</span></div>

<div v-click="1" class="mot-cfd">
<div class="name">one resolved condition costs</div>
<div class="row"><span class="big2">≈ 4,000</span> CPU-h &nbsp;<span class="cap">(our L10)</span></div>
<div class="quote">"highly time consuming and not practically feasible for a large number of operating conditions"<span class="who"> — Zhan 2019</span></div>
</div>
</div>

<div class="mot-right" v-click="2">
<img src="/lit_gap.png" />
</div>

<style>
.mot-left { position: absolute; left: 34px; top: 24px; width: 520px; }
.mot-right { position: absolute; right: 22px; top: 40px; width: 420px; }
.mot-right img { width: 420px; }
.mot-cost { display: flex; align-items: center; gap: 22px; margin: 10px 0 0 4px; }
.mot-cost .big { font-size: 46px; font-weight: 300; color: var(--ink); }
.mot-cost .u { font-size: 20px; color: var(--mute); }
.mot-cost .arrow { font-size: 34px; color: var(--mute); }
.mot-eq .katex-display { margin: 10px 0 4px; }
.mot-eq .katex { font-size: 1.0em; }
.quote { font-size: 16px; font-style: italic; color: var(--ink); margin-top: 4px; }
.quote .who { font-style: normal; color: var(--mute); }
.mot-cfd { margin-top: 18px; }
.mot-cfd .row { font-size: 22px; color: var(--ink); margin-top: 2px; }
.mot-cfd .big2 { font-size: 36px; font-weight: 300; }
.cap { color: var(--mute); font-size: 15px; line-height: 1.2; }
</style>

<!--
Cultivated meat is still far from the price of beef. Negulescu and co-workers estimate 25 to 35 dollars per kilo in stirred tanks of 42 to 211 cubic metres. To compete with beef, production must drop below 9 dollars per kilo. Media dominate the cost, but the bioreactor sets the cell density and the scale.

The bioreactor must do three things at once. The shear stress must stay below the damage range, 1 to 100 pascal for mammalian cells. The oxygen transfer must supply the oxygen uptake of the cells. The mixing time must be shorter than one over kLa, so the oxygen spreads before the cells use it. These pull against each other: more agitation helps oxygen and mixing, but raises shear.

[click] Simulations can map this trade-off before experiments. But one resolved condition costs about 4,000 CPU hours at our finest grid. Zhan and co-workers say this directly for their wave bag.

[click] What earlier work did. Zhan, Seidel and Kim simulate rocking bags with a resolved interface, each on one grid, with no error bars. Cantarero-Rivera corrects a coarse mesh with a neural network, in a spinner flask, again with no error bars, and it does not generalise across speeds. Kiviat computes all three quantities, but in stirred tanks with RANS, which they say under-predicts shear. Martens uses a multi-fidelity Gaussian process, but the fidelities are culture scales from an ODE model, not grid levels. We do all five: rocking bag, resolved two-phase flow, all three quantities, grid levels used as data, and error bars.
-->

---

# Introduction

---
clicks: 3
---

<div class="g4">
<div class="cell">
<div class="head"><span class="name">shear stress</span><div class="eq4">

$\tau = \mu\left(\dfrac{\partial u}{\partial y} + \dfrac{\partial v}{\partial x}\right)$

</div></div>
<div class="vid">
  <SwapVideo loop="/tau_loop.mp4" rel="/oxy_rel.mp4" :on="$clicks >= 1" />
  <div class="cbar" :class="{ shown: $clicks >= 1 }"><span class="bar"></span><span class="eq"><i>c</i></span></div>
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
  <div class="cbar" :class="{ shown: $clicks >= 3 }"><span class="bar"></span><span class="eq"><i>c</i></span></div>
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
.g4 { position: absolute; inset: 26px 34px 20px 34px; display: grid; grid-template-columns: 1fr 1fr; grid-template-rows: 1fr 1fr; gap: 14px 40px; }
.cell { position: relative; display: flex; flex-direction: column; min-height: 0; }
.head { display: flex; align-items: center; gap: 18px; height: 62px; }
.head .name { font-size: 26px; }
.eq4 p, .note4 p { margin: 0; }
.eq4 .katex { font-size: 1.32em; color: var(--ink); }
.vid { position: relative; flex: 1; display: flex; flex-direction: column; justify-content: center; }
.plot { height: 178px; margin: 0 !important; align-self: flex-start; }
.note4 { position: absolute; left: 64px; top: 66px; }
.note4 .katex { font-size: 0.95em; color: var(--mute); }
.cbar { display: flex; align-items: center; gap: 8px; align-self: flex-end; margin-top: 8px; opacity: 0; transition: opacity 0.8s; }
.cbar.shown { opacity: 1; }
.cbar .bar { width: 110px; height: 8px; border-radius: 4px; background: linear-gradient(90deg, #eef3f8, #b9dcee, #5fb7d6, #1f6fc0, #1d2a78); }
.cbar .eq { font-size: 20px; color: var(--ink); }
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
clicks: 3
---

<div class="mdl-left">
<div class="name">two-phase Navier–Stokes, in the frame of the bag</div>

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

<div class="mdl-grid">
<img src="/quadtree.png" />
<div>

$${\color{#0072B2}h_8 = \tfrac{L}{2^{8}}} \quad {\color{#009E73}h_9 = \tfrac{L}{2^{9}}} \quad {\color{#CC79A7}h_{10} = \tfrac{L}{2^{10}}}$$

<div class="cap">Basilisk, adaptive quadtree</div>
</div>
</div>

</div>

<div class="mdl-right" v-click="1">
<div class="name">multi-fidelity</div>

$${\color{#CC79A7}f_H(x)} = \rho\,{\color{#0072B2}f_L(x)} + b + \delta(x), \qquad \delta \sim \mathcal{GP}(0,\, s^2 k_\ell) $$

<div class="mdl-fig">
  <img v-if="$clicks <= 1" src="/mf_a.png" />
  <img v-else-if="$clicks == 2" src="/mf_b.png" />
  <video v-else src="/mf_iter.mp4" autoplay muted playsinline />
</div>
<div class="mdl-cost"><span style="color:#CC79A7">one L10 run ≈ <b>4,000</b> CPU-h</span> &nbsp;·&nbsp; <span style="color:#0072B2">one L8 run ≈ <b>60</b> CPU-h</span></div>

<div v-click="3">

$$\theta^\star = \arg\min_\theta\, \left[-\log p(\mathbf{y}_H \mid \theta)\right], \qquad \theta = (\rho,\, b,\, \ell,\, s,\, \sigma_n)$$

</div>
</div>

<style>
.mdl-left { position: absolute; left: 34px; top: 24px; width: 520px; }
.mdl-right { position: absolute; right: 26px; top: 24px; width: 400px; }
.mdl-left .katex-display, .mdl-right .katex-display { margin: 6px 0; }
.mdl-left .katex-display .katex { font-size: 0.95em; }
.mdl-row { display: flex; justify-content: space-around; margin-top: 4px; text-align: center; }
.mdl-grid { display: flex; align-items: center; gap: 14px; margin-top: 8px; }
.mdl-grid img { height: 92px; }
.mdl-cost { margin-top: 2px; font-size: 15px; white-space: nowrap; color: var(--ink); }
.mdl-fig img, .mdl-fig video { width: 400px; display: block; }
.mdl-right .katex-display .katex { font-size: 0.76em; }
.mdl-right .katex-display { overflow: visible; }
.cap { color: var(--mute); font-size: 18px; white-space: nowrap; }
</style>

<!--
The solver is Basilisk. We solve the two-phase Navier–Stokes equations in the frame of the bag, so the walls do not move. The price is the frame terms, in orange: the angular acceleration, Coriolis and centrifugal terms. The volume fraction alpha marks the liquid. At the interface the oxygen follows Henry's law. The mesh is a quadtree. L8, L9 and L10 mean a cell size of L over 2 to the 8, 9 and 10.

[click] Multi-fidelity. Blue dots: many cheap L8 runs. Pink dots: few expensive L10 runs. One L10 run costs about 4,000 CPU hours (3,600 measured at 32.5 rpm; more at low rpm). One L8 run costs about 60. So we can afford few L10 runs.

[click] A kernel ridge regression fits the cheap data: f_L.

[click] The expensive data are rho times f_L, plus b, plus a Gaussian process delta. We fit rho, b and the GP parameters by minimising the negative log-likelihood. Each pink curve is one step of the optimiser. The newest step is the strong line. On the right, the bold parameters are the ones that change at that step. At the end we get the prediction and its 95 % band. The data here are synthetic. Our data are in the results. We chain this: L8 to L9, then L9 to L10.
-->

---

# Results 1

---

# Results 2

---

# Conclusions

---

# Questions?
