"""Rewrite the Methods slide block of slides.md (the '# Methods' placeholder or a previous version) from the template below."""
from pathlib import Path

NEW = '''<div class="absolute left-10 top-8" style="width: 400px">
  <div v-click>
    <div class="flex items-baseline justify-between"><span class="name">Basilisk</span><span class="cap">Navier–Stokes + VOF</span></div>
    <img src="/met_bag.png" style="height: 190px" />
  </div>
</div>

<div class="absolute right-10 top-8" style="width: 470px">
  <div v-click>
    <div class="flex items-baseline justify-between"><span class="name">release</span><span class="cap">∂c/∂t + u·∇c = ∇·(D∇c)</span></div>
    <img src="/met_timeline.png" style="height: 190px" />
  </div>
</div>

<div class="absolute left-10 right-10" style="top: 320px">
  <div v-click>
    <div class="flex items-baseline justify-between"><span class="name">multi-fidelity</span><span class="cap">log QoI : L8 → L9 → L10</span></div>
    <img src="/met_chain.png" style="height: 175px" />
  </div>
</div>

<!--
[click] The solver is Basilisk. It solves the Navier–Stokes equations with a volume-of-fluid method. The bag is an embedded boundary. We rock it by plus and minus seven degrees. Oxygen enters through the free surface.

[click] We settle the flow for 80 rocking cycles. Then we release the tracer in the top half of the liquid and start the oxygen. Both follow an advection–diffusion equation. We measure the curves after the release.

[click] The method is multi-fidelity. L8 is cheap and has many points. L9 has fewer points. L10 is expensive and has the fewest. A Gaussian process chain on the log of each quantity carries the information from L8 to L9 to L10.
-->
'''
p = Path(__file__).resolve().parent.parent / "slides.md"
s = p.read_text()
if "# Methods" in s:
    a = s.index("# Methods"); b = s.index("---", a)
    p.write_text(s[:a] + NEW + "\n" + s[b:])
else:
    raise SystemExit("Methods placeholder not found")
print("ok")
