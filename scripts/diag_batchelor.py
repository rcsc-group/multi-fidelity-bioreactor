"""Kolmogorov and Batchelor scales from our own dissipation data, compared with the cell size of each level.

eta_K = (nu^3 / eps)^(1/4); eta_B = eta_K / sqrt(Sc) = (nu D^2 / eps)^(1/4), Sc = nu / D (Batchelor 1959; the
standard definition). eps (W/kg) = edr_mean_t / rho, with edr_mean_t from experiments/multifidelity/
hydro_dataset_v2.json: the DOMAIN-mean dissipation (W/m^3, liquid + air) over the last 3 cycles before release.
Because the air carries little dissipation, the liquid-mean eps is larger by up to 1 / (liquid area fraction); the
scales move as eps^(-1/4), so this uncertainty is a factor <= (1/f_liq)^(1/4) (printed).
Cell size dx = L_bio / 2^L with L_bio = 0.25 m (src/BioReactor.c: L0 = 1 in units of L_bio).

Usage: uv run python scripts/diag_batchelor.py
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NU, RHO, L_BIO = 1.0e-6, 1.0e3, 0.25
D = {"tracer": 0.44e-9, "oxygen": 1.9e-9}  # src/BioReactor.c D_tracer_1, D_oxy_1
F_LIQ = 0.286 / 1.0  # liquid area / domain area is NOT known here; Kim's Ly = 0.286 is a domain height, so this is
# only a placeholder bound used for the printed uncertainty factor; see docstring.


def main():
    rows = json.load(open(ROOT / "experiments/multifidelity/hydro_dataset_v2.json"))
    print(f"{'L':>3} {'rpm':>5} {'eps W/kg':>9} {'eta_K um':>9} {'eta_B tr um':>11} {'eta_B O2 um':>11} "
          f"{'dx um':>7} {'dx/eta_K':>8} {'dx/eta_B tr':>11}")
    out = []
    for r in sorted(rows, key=lambda r: (r["level"], r["rpm"])):
        if r["rpm"] not in (17.5, 25.0, 32.5, 37.5):
            continue
        eps = r["edr_mean_t"] / RHO
        eta_k = (NU ** 3 / eps) ** 0.25
        eta_b = {k: eta_k / (NU / d) ** 0.5 for k, d in D.items()}
        dx = L_BIO / 2 ** r["level"]
        out.append(dict(level=r["level"], rpm=r["rpm"], run=r["run"], eps=eps, eta_K=eta_k, eta_B=eta_b, dx=dx))
        print(f"{r['level']:>3} {r['rpm']:>5g} {eps:9.2e} {eta_k*1e6:9.0f} {eta_b['tracer']*1e6:11.1f} "
              f"{eta_b['oxygen']*1e6:11.1f} {dx*1e6:7.0f} {dx/eta_k:8.2f} {dx/eta_b['tracer']:11.0f}")
    print(f"uncertainty from domain- vs liquid-mean eps: scales smaller by up to {(1/F_LIQ)**0.25:.2f}x (placeholder)")
    json.dump(out, open(ROOT / "experiments/multifidelity/diag_batchelor.json", "w"), indent=1)


if __name__ == "__main__":
    main()
