"""Step 4 (lead, 2026-11-05): LOO residuals of a stationary 1-D GP in rpm along theta 4 (L8 runs, rpm 20-27.5, read from
diag_tau_l8_local output). Fixed length scales scanned (no hyperparameter fit), noise = 0.017 ln (tau floor).
Prints LOO residuals (ln) near the feature (22-24 rpm) and far from it (20-21, 26-27.5)."""
import numpy as np

raw = """20.000 0.00680
20.625 0.00701
21.250 0.00710
21.750 0.00715
21.875 0.00712
22.000 0.00734
22.250 0.00734
22.500 0.00723
22.750 0.00673
23.000 0.00692
23.125 0.00712
23.250 0.00726
23.750 0.00746
24.375 0.00746
25.000 0.00745
25.625 0.00775
26.250 0.00782
26.875 0.00791
27.500 0.00790"""
d = np.array([[float(a), float(b)] for a, b in (ln.split() for ln in raw.splitlines())])
x, y = d[:, 0], np.log(d[:, 1])
y = y - y.mean()
noise = 0.017
print(f"n={len(x)}, sd(ln) of data={y.std(ddof=1):.4f}, noise={noise}")
for ell in (0.5, 1.0, 2.0, 4.0):
    K = np.exp(-0.5 * ((x[:, None] - x[None, :]) / ell) ** 2) * (y.var() + 1e-12) + noise ** 2 * np.eye(len(x))
    Ki = np.linalg.inv(K)
    loo = (Ki @ y) / np.diag(Ki)            # y_i - mu_{-i}
    near = (x >= 22.0) & (x <= 24.0)
    far = (x < 21.0) | (x > 26.0)
    print(f"ell={ell:4.1f} rpm: LOO rms near 22-24 = {np.sqrt(np.mean(loo[near] ** 2)):.4f}  far = {np.sqrt(np.mean(loo[far] ** 2)):.4f}  "
          f"max|LOO| near = {np.max(np.abs(loo[near])):.4f} at {x[near][np.argmax(np.abs(loo[near]))]}")
