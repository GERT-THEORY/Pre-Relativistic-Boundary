"""
GERT Paper 3 — Emergence boundary with Saha equilibrium
=========================================================
Section 3.5–3.6

Computes the metric emergence parameter Ξ(α) = λ_γ / d_ph using the
Saha equilibrium ionization fraction, and determines α_em from Ξ = 1.

Result: α_em = −3.047 ± 0.001 (z ≈ 1114, T ≈ 3037 K)
"""
import numpy as np
from scipy.optimize import brentq
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from gert_utils_p3 import *

print("=" * 72)
print("  GERT Paper 3 — Ξ(α) with Saha Equilibrium")
print("=" * 72)

# ── Analytical estimates ────────────────────────────────────────────────
A_opacity = c * n_b0 * sigma_T / (H0 * np.sqrt(Omega_r0))

print(f"\n  Fundamental scales:")
print(f"    n_b0        = {n_b0:.4e} m⁻³")
print(f"    λ_0         = {lambda_0:.4e} m = {lambda_0/Mpc:.0f} Mpc")
print(f"    Ω_r,0       = {Omega_r0:.4e}")
print(f"    A (opacity–Hubble scale) = {A_opacity:.6f}")

a_rec = 1/1101
alpha_rec = z_to_alpha(1100)
lam_rec = lambda_0 * a_rec**3
d_ph_approx = c * a_rec**2 / (H0 * np.sqrt(Omega_r0))
Xi_before = lam_rec / d_ph_approx

print(f"\n  At α = {alpha_rec:.4f} (z = 1100, CMB):")
print(f"    λ_γ(X_e=1)  = {lam_rec:.3e} m = {lam_rec/Mpc*1e3:.2f} kpc")
print(f"    d_ph (rad)   = {d_ph_approx:.3e} m = {d_ph_approx/Mpc*1e3:.1f} kpc")
print(f"    Ξ(before)    = {Xi_before:.4e}  (≪ 1: NOT readable)")
print(f"    Ξ(after)     ≈ {Xi_before*1e4:.1f}  (≫ 1: readable)")

# ── Ξ(α) profile ───────────────────────────────────────────────────────
print(f"\n{'='*72}")
print(f"  Ξ(α) profile across recombination")
print(f"{'='*72}")
print(f"\n  {'α':<10} {'x':<12} {'z':<8} {'T(K)':<8} "
      f"{'X_e':<12} {'λ_γ(kpc)':<10} {'d_ph(kpc)':<10} {'Ξ':<12}")
print(f"  {'-'*90}")

alpha_scan = [-4.0, -3.5, -3.3, -3.2, -3.15, -3.10, -3.08,
              -3.06, -3.05, -3.04, -3.03, -3.02, -3.00,
              -2.95, -2.90, -2.80, -2.50, -2.00, 0.00]

for alpha in alpha_scan:
    z = alpha_to_z(alpha)
    x = alpha_to_x(alpha)
    T = alpha_to_T(alpha)
    Xe = Xe_saha(alpha)
    lam = lambda_gamma(alpha) / Mpc * 1e3   # kpc
    dph = d_ph(alpha) / Mpc * 1e3            # kpc
    xi = Xi(alpha)

    mark = "  ← Ξ ≈ 1" if 0.3 < xi < 3 else ""
    print(f"  {alpha:<10.4f} {x:<12.2f} {z:<8.0f} {T:<8.0f} "
          f"{Xe:<12.4e} {lam:<10.2f} {dph:<10.1f} {xi:<12.4e}{mark}")

# ── Find α_em ──────────────────────────────────────────────────────────
print(f"\n{'='*72}")
print(f"  Emergence boundary: Ξ(α_em) = 1")
print(f"{'='*72}")

alpha_em = brentq(lambda a: Xi(a) - 1.0, -3.3, -2.5, xtol=1e-12)
z_em  = alpha_to_z(alpha_em)
x_em  = alpha_to_x(alpha_em)
T_em  = alpha_to_T(alpha_em)
Xe_em = Xe_saha(alpha_em)

print(f"\n  α_em  = {alpha_em:.6f}")
print(f"  x_em  = {x_em:.4f}  (log₁₀ρ at emergence)")
print(f"  z_em  = {z_em:.1f}")
print(f"  T_em  = {T_em:.0f} K  ({T_em*k_B/eV:.4f} eV)")
print(f"  X_e   = {Xe_em:.6f}")

# Verify
xi_check = Xi(alpha_em)
print(f"\n  Verification: Ξ(α_em) = {xi_check:.8f}  ✓")

# Domain
print(f"\n  GERT relativistic domain (Saha):")
print(f"    α_em   = {alpha_em:.2f}")
print(f"    α_crit = {alpha_crit:.2f}")
print(f"    Span   = {alpha_crit - alpha_em:.1f} decades in α")
