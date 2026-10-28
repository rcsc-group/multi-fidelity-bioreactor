"""Write fonts_test.md: the Definitions slide in light mode, once per candidate font class (fA, fB, fC defined in style.css)."""
from pathlib import Path

CELL = '''
<div class="absolute left-10 top-8" style="width: 440px">
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
    <img src="/def_tau_hist_light.png" class="absolute" style="right: 0; top: 36px; height: 86px" />
  </div>
  <div v-click style="margin-top: 60px">
    <div class="flex items-baseline justify-between"><span class="name">dissipation rate</span><span class="cap">ε = 2 μ S : S</span></div>
    <img src="/def_edr_light.png" style="width: 470px" />
  </div>
</div>
'''
head = "---\ntheme: default\naspectRatio: 16/9\ncanvasWidth: 1024\nmdc: true\ncolorSchema: light\ntransition: fade\n---\n"
parts = [head]
for i, cls in enumerate(("fA", "fB", "fC")):
    sep = "" if i == 0 else "\n---\n"
    parts.append(f'{sep}\n<div class="{cls} absolute inset-0">\n{CELL}\n</div>\n')
Path(__file__).resolve().parent.parent.joinpath("fonts_test.md").write_text("\n".join(parts))
print("ok")
