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
---

# Motivation

---

# Introduction

---

<div class="absolute left-10 top-8" style="width: 440px">
  <div v-click>
    <div class="flex items-baseline justify-between"><span class="name">kLa</span><span class="cap eq">dC*/dt = kLa (1 − C*)</span></div>
    <img src="/def_kla_light.png" style="height: 205px" />
  </div>
  <div v-click style="margin-top: 18px">
    <div class="flex items-baseline justify-between"><span class="name">mixing time</span><span class="cap eq">χ = 1 − σ² / σ²max</span></div>
    <img src="/def_mix_light.png" style="height: 200px" />
  </div>
</div>

<div class="absolute right-10 top-8" style="width: 470px">
  <div v-click class="relative">
    <div class="flex items-baseline justify-between"><span class="name">shear stress</span><span class="cap eq">τ = μ (∂u/∂y + ∂v/∂x)</span></div>
    <img src="/def_tau_light.png" style="width: 470px" />
    <img src="/def_tau_hist_light.png" style="height: 78px; margin: 6px 0 0 auto" />
  </div>
  <div v-click style="margin-top: 26px">
    <div class="flex items-baseline justify-between"><span class="name">dissipation rate</span><span class="cap eq">ε = 2 μ S : S</span></div>
    <img src="/def_edr_light.png" style="width: 470px" />
  </div>
</div>

<!--
[click] kLa is the speed of oxygen transfer. Oxygen goes in at the surface. We fit the rise of the saturation over one rocking period. We read the value at ten percent saturation.

[click] Mixing time is the time to mix a tracer. At the start the tracer fills the top half. The variance of the tracer falls. We stop the clock at fifty percent mixing.

[click] Shear stress hurts the cells. We take the 95th percentile of the stress over the liquid.

[click] The dissipation rate is the energy loss in the liquid. It is the second stress measure.

All four come from one run of 32.5 rpm and 7 degrees. The values are kLa 47 per hour and mixing time 11.2 seconds at the lowest fidelity.
-->

---

<div class="absolute left-10 top-8" style="width: 400px">
  <div v-click>
    <div class="flex items-baseline justify-between"><span class="name">Basilisk</span><span class="cap">Navier–Stokes + VOF</span></div>
    <img src="/met_bag.png" style="height: 150px; margin: 6px auto 0" />
  </div>
  <div v-click style="margin-top: 26px">
    <div class="flex items-baseline justify-between"><span class="name">release</span><span class="cap eq">∂c/∂t + u·∇c = ∇·(D∇c)</span></div>
    <img src="/met_timeline.png" style="height: 150px; margin-top: 6px" />
  </div>
</div>

<div class="absolute right-10 top-8" style="width: 500px">
  <div v-click>
    <div class="name">multi-fidelity</div>
    <img src="/mf_schematic.png" style="width: 500px" />
    <div class="eq" style="margin-top: 8px; text-align: center; font-size: 32px; white-space: nowrap">
      <span style="color: #1d2a78"><i>f</i><sub>H</sub>(<i>x</i>)</span> = <span style="color: #7a8797"><i>ρ f</i><sub>L</sub>(<i>x</i>)</span> + <span style="color: #e8590c"><i>δ</i>(<i>x</i>)</span>
    </div>
    <div class="cap eq" style="text-align: center; margin-top: 2px"><i>δ</i> ∼ GP</div>
  </div>
</div>

<!--
[click] The solver is Basilisk. It solves the Navier–Stokes equations with a volume-of-fluid method. The bag is an embedded boundary. We rock it by plus and minus seven degrees. Oxygen enters through the free surface.

[click] We settle the flow for 80 rocking cycles. Then we release the tracer in the top half of the liquid and start the oxygen. Both follow an advection–diffusion equation. We measure the curves after the release.

[click] The method is multi-fidelity. It follows Yi and co-workers. The light dots are cheap runs. There are many. They give the dashed low-fidelity surrogate. The grey curve is that surrogate after a linear map. The dark dots are expensive runs. There are few. The orange sticks are the gap between the grey curve and the dark dots. A Gaussian process models this gap. The dark curve is the sum. We apply the step from L8 to L9, and again from L9 to L10. We work on the log of each quantity. The figure is a sketch, not our data. Our data are in the results.
-->

---

# Results 1

---

# Results 2

---

# Conclusions

---

# Questions?
