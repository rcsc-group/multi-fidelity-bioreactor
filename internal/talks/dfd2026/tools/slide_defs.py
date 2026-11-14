"""Rewrite the Definitions slide block of slides.md from the template below (light mode, text-free figures in ../public/def_*_light.png)."""
from pathlib import Path

NEW = '''<div class="absolute left-10 top-8" style="width: 440px">
  <div v-click>
    <div class="flex items-baseline justify-between"><span class="name">kLa</span><span class="cap">dC*/dt = kLa (1 − C*)</span></div>
    <img src="/def_kla_light.png" style="height: 205px" />
  </div>
  <div v-click style="margin-top: 18px">
    <div class="flex items-baseline justify-between"><span class="name">mixing time</span><span class="cap">χ = 1 − σ² / σ²max</span></div>
    <img src="/def_mix_light.png" style="height: 200px" />
  </div>
</div>

<div class="absolute right-10 top-8" style="width: 470px">
  <div v-click class="relative">
    <div class="flex items-baseline justify-between"><span class="name">shear stress</span><span class="cap">τ = μ (∂u/∂y + ∂v/∂x)</span></div>
    <img src="/def_tau_light.png" style="width: 470px" />
    <img src="/def_tau_hist_light.png" style="height: 78px; margin: 6px 0 0 auto" />
  </div>
  <div v-click style="margin-top: 26px">
    <div class="flex items-baseline justify-between"><span class="name">dissipation rate</span><span class="cap">ε = 2 μ S : S</span></div>
    <img src="/def_edr_light.png" style="width: 470px" />
  </div>
</div>
'''
p = Path(__file__).resolve().parent.parent / "slides.md"
s = p.read_text()
a = s.index('<div class="absolute left-10 top-6 cap"') if '<div class="absolute left-10 top-6 cap"' in s else s.index('<div class="absolute left-10 top-8"')
b = s.index("<!--", a)
p.write_text(s[:a] + NEW + "\n" + s[b:])
print("ok")
