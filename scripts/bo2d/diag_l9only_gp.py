"""Why is the L9-only kLa GP flat with low latent sd (2026-11-09)? Print the posterior over (length, noise ratio),
the implied noise sd in ln units, and the data spread. Hypothesis: posterior mass at the largest noise ratio / long
length -> data treated as noise around a constant."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import bo, loop as L  # noqa: E402
K = L.arm_data("l9only", L.QK)[9]
y = np.log(K.y.to_numpy(float))
m = bo.fit_single(K, L.QK)
print("n", len(y), "ln y sd", y.std().round(3), "range", np.exp(y.min()).round(1), np.exp(y.max()).round(1))
print("noise floor (ln)", bo.C.noise_floor(L.QK))
print("grid: length", np.unique(m.l).round(3), "eta", np.unique(m.eta).round(3))
w = m.w / m.w.sum()
for name, v in (("length", m.l), ("eta", m.eta)):
    print(f"posterior mean {name} {np.sum(w * v):.3f}; mass at max {w[v == v.max()].sum():.2f}; at min {w[v == v.min()].sum():.2f}")
mu, V, nz = bo.latent(m, K[["rpm", "theta"]].to_numpy(float))
print("noise sd (ln)", np.sqrt(nz).round(3), "latent sd range", np.sqrt(V).min().round(3), np.sqrt(V).max().round(3))
print("residual sd (ln)", np.std(y - mu).round(3))
