"""Rewrite the Methods block of slides.md: Basilisk and release (small, left), ONE multi-fidelity hero figure (right). Light mode."""
from pathlib import Path

NEW = '''<div class="absolute left-10 top-8" style="width: 400px">
  <div v-click>
    <div class="flex items-baseline justify-between"><span class="name">Basilisk</span><span class="cap">Navier–Stokes + VOF</span></div>
    <img src="/met_bag.png" style="height: 150px; margin: 6px auto 0" />
  </div>
  <div v-click style="margin-top: 26px">
    <div class="flex items-baseline justify-between"><span class="name">release</span><span class="cap">∂c/∂t + u·∇c = ∇·(D∇c)</span></div>
    <img src="/met_timeline.png" style="height: 150px; margin-top: 6px" />
  </div>
</div>

<div class="absolute right-10 top-8" style="width: 500px">
  <div v-click>
    <div class="name">multi-fidelity</div>
    <img src="/mf_hero_v3.png" style="width: 500px" />
    <div class="cap" style="margin-top: 4px; display: flex; gap: 16px; justify-content: center">
      <span style="color: #9fc3e8">● L8</span><span style="color: #4f93d6">● L9</span><span style="color: #1d2a78">● L10</span><span>▮ 95 % band</span>
    </div>
  </div>
</div>

'''
p = Path(__file__).resolve().parent.parent / "slides.md"
s = p.read_text()
k = s.index('<span class="name">Basilisk</span>')
a = s.rindex('<div class="absolute left-10 top-8"', 0, k)
b = s.index("<!--", k)
p.write_text(s[:a] + NEW + s[b:])
print("ok")
