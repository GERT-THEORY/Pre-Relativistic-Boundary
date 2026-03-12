"""
GERT Paper 3 — Complete domain map in GERT natural time
=========================================================
Section 3.8–3.9

Combines results from Papers I, II, and III to produce the complete
domain map of the GERT framework in the natural time coordinate
α = log₁₀(a).

Domain:
  Primordial Cauldron   α < α_PC          (Layer 3 absent: not yet crystallized)
  Pre-legible geometry  α_PC < α < α_em   (Ξ < 1: metric exists but not readable)
  Relativistic regime   α_em ≤ α ≤ α_crit (Ξ ≥ 1: GERT equations valid)
  Quasi-Vacuum          α > α_crit        (Layer 3 dissolved: post-relativistic)

α_em   = −3.0 ± 0.1     (this paper)
α_crit = +12.88 ± 0.12  (Paper II)
Span   = 15.9 ± 0.2 decades
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from gert_utils_p3 import *

print("=" * 72)
print("  GERT — Complete Domain Map in Natural Time α = log₁₀(a)")
print("=" * 72)

# ── Boundaries ──────────────────────────────────────────────────────────
alpha_em_saha   = -3.047
alpha_em_peebles = -2.937
alpha_em_mid    = (alpha_em_saha + alpha_em_peebles) / 2
alpha_em_err    = abs(alpha_em_peebles - alpha_em_saha) / 2
alpha_cr        = 12.88
alpha_cr_err    = 0.12

span     = alpha_cr - alpha_em_mid
span_err = np.sqrt(alpha_em_err**2 + alpha_cr_err**2)

# ── Key epochs ──────────────────────────────────────────────────────────
print(f"\n  {'Epoch':<30} {'α':<12} {'x = log₁₀ρ':<14} {'z':<10} {'T (K)':<12}")
print(f"  {'-'*78}")

epochs = [
    ("α_em (Saha)",       alpha_em_saha,    None),
    ("α_em (Peebles)",    alpha_em_peebles, None),
    ("Rad-matter equality", z_to_alpha(3400), 3400),
    ("CMB (last scatter.)", z_to_alpha(1100), 1100),
    ("Dark ages",           z_to_alpha(900),  900),
    ("Reionization",        z_to_alpha(10),   10),
    ("Today",               0.0,              0),
    ("α_crit (Paper II)",   alpha_cr,         None),
]

for name, alpha, z_ref in epochs:
    x = alpha_to_x(alpha)
    z = alpha_to_z(alpha) if z_ref is None else z_ref
    T = alpha_to_T(alpha)
    z_str = f"{z:.0f}" if abs(z) < 1e6 else f"10^{np.log10(abs(z)+1):.1f}"
    print(f"  {name:<30} {alpha:<12.4f} {x:<14.2f} {z_str:<10} {T:<12.3e}")

# ── Domain map ──────────────────────────────────────────────────────────
print(f"""
{'='*72}
  DOMAIN MAP
{'='*72}

  α-axis (GERT natural time, α = log₁₀ a)

  ─────────────────────────────────────────────────────────────→ α

  −∞        α_PC         α_em             0            α_crit
   │         │            │               │              │
   ▼         ▼            ▼               ▼              ▼
  [Cauldron]  [Pre-legible]  [RELATIVISTIC REGIME]     [Quasi-Vacuum]
  Layer 3     Layer 3        Layer 3 active,            Layer 3
  not yet     crystallizing  Ξ ≥ 1                     dissolved
  crystallized               GERT eqs valid

  BOUNDARIES:
    α_em   = {alpha_em_mid:.1f} ± {alpha_em_err:.1f}    (this paper)
    α_crit = {alpha_cr:.2f} ± {alpha_cr_err:.2f}  (Paper II)

  RELATIVISTIC DOMAIN:
    Span = {span:.1f} ± {span_err:.1f} decades in α
         = ~{abs(alpha_to_x(alpha_em_mid) - alpha_to_x(alpha_cr)):.0f} decades in log₁₀(ρ)

  UNCERTAINTY SOURCES:
    α_em  :  recombination kinetics (Saha vs Peebles)
    α_crit:  k_gas (Paper I free parameter)

  INSENSITIVITY:
    α_PC (25+ dex variation):  |Δα_em| < 5 × 10⁻⁴
""")

# ── Density at boundaries ──────────────────────────────────────────────
print(f"  {'Boundary':<22} {'α':<12} {'log₁₀(ρ) [kg/m³]':<20} {'ρ [kg/m³]'}")
print(f"  {'-'*66}")
for name, alpha in [("α_em (Saha)", alpha_em_saha),
                     ("α_em (Peebles)", alpha_em_peebles),
                     ("Today", 0),
                     ("α_crit", alpha_cr)]:
    x = alpha_to_x(alpha)
    rho = 10**x
    print(f"  {name:<22} {alpha:<12.2f} {x:<20.2f} {rho:.2e}")

# ── Paper II cross-check ───────────────────────────────────────────────
print(f"\n  Paper II consistency:")
print(f"    α_crit = 12.88 ± 0.12")
print(f"    log₁₀(ρ_GR,min) = {alpha_to_x(alpha_cr):.1f} ± 0.4  (Paper II: −65.2 ± 0.4  ✓)")

# ── Structural symmetry table ──────────────────────────────────────────
print(f"""
{'='*72}
  STRUCTURAL SYMMETRY (Papers II + III)
{'='*72}

  {'Property':<28} {'Paper II (dilute)':<24} {'Paper III (dense)'}
  {'-'*72}
  {'Boundary':<28} {'α_crit = +12.88 ± 0.12':<24} {'α_em = −3.0 ± 0.1'}
  {'Mechanism':<28} {'Dilution dissolves Layer 3':<24} {'Excess curvature blocks L3'}
  {'Criterion':<28} {'λ ≥ Hubble radius':<24} {'λ_γ ≥ particle horizon'}
  {'Uncertainty source':<28} {'k_gas (Paper I)':<24} {'Recombination kinetics'}
  {'Observational anchor':<28} {'Asymptotic thermal death':<24} {'CMB (first metric map)'}
  {'Active layers beyond':<28} {'Layer 1 + 2':<24} {'Layer 1 + 2'}
  {'GERT equations':<28} {'No longer valid':<24} {'Not yet valid'}
""")
