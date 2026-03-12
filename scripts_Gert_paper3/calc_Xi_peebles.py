"""
GERT Paper 3 — Emergence boundary with Peebles recombination
==============================================================
Section 3.6

Solves the Peebles three-level atom equation for X_e(z) and recomputes
the metric emergence parameter Ξ(α). The Peebles C-factor captures the
kinetic bottleneck of recombination: most atoms reaching n=2 are
re-photoionized before decaying to the ground state, slowing recombination
and raising the freeze-out residual to X_e ~ 3×10⁻⁴.

Result: α_em = −2.937 ± 0.001 (z ≈ 865, T ≈ 2358 K)

Combined with Saha (α_em = −3.047): α_em = −3.0 ± 0.1

Reference: Peebles (1968), ApJ 153, 1;
           Seager, Sasselov & Scott (1999), ApJ 523, L1
"""
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import brentq
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from gert_utils_p3 import *

# ══════════════════════════════════════════════════════════════════════════
#  Peebles three-level atom ODE
# ══════════════════════════════════════════════════════════════════════════

def C_peebles(Xe, z):
    """Peebles C-factor: fraction of atoms in n=2 that reach ground state
    before being re-photoionized.

    C = (Λ_2s + Λ_Lyα) / (Λ_2s + Λ_Lyα + β_B)

    where Λ_2s is the 2s→1s two-photon rate, Λ_Lyα is the net
    Lyman-α escape rate, and β_B is the photoionization rate from n=2.
    """
    T   = T0 * (1 + z)
    kT  = k_B * T
    nH  = n_H0 * (1 + z)**3
    n1s = nH * (1 - Xe)

    # Photoionization rate from n=2 (detailed balance with α_B)
    thermal = (m_e * kT / (2 * np.pi * hbar**2))**1.5
    beta_B  = alpha_B(T) * thermal * np.exp(-B_2 / kT)

    Hz = H_z(z)

    if n1s < 1.0:
        return 1.0

    # Net Lyman-α escape rate (Sobolev approximation)
    Lam_Lya = Hz * 8 * np.pi / (n1s * lambda_Lya**3)

    return (Lambda_2s + Lam_Lya) / (Lambda_2s + Lam_Lya + beta_B)


def dXe_dz_peebles(Xe, z):
    """Peebles ODE: dX_e/dz.

    dX_e/dt = C · [β_B · exp(−E_Lyα/kT) · (1−X_e) − α_B · n_H · X_e²]

    The factor exp(−E_Lyα/kT) ensures detailed balance: combined with
    β_B ∝ exp(−B_2/kT), it gives the correct ground-state ionization
    rate ∝ exp(−B_1/kT) since B_2 + E_Lyα = B_1.

    Using dt/dz = −1/[(1+z)H]:
    dX_e/dz = C/[(1+z)H] · [α_B · n_H · X_e² − β_B · exp(−E_Lyα/kT) · (1−X_e)]
    """
    Xe = np.clip(Xe, 1e-8, 1.0 - 1e-10)

    T   = T0 * (1 + z)
    kT  = k_B * T
    nH  = n_H0 * (1 + z)**3
    Hz  = H_z(z)
    aB  = alpha_B(T)

    # Rates
    thermal = (m_e * kT / (2 * np.pi * hbar**2))**1.5
    beta_B  = aB * thermal * np.exp(-B_2 / kT)

    ion_rate = beta_B * np.exp(-E_Lya / kT) * (1 - Xe)
    rec_rate = aB * nH * Xe**2

    Cp = C_peebles(Xe, z)

    # dXe/dz > 0 when recombining (z decreasing → Xe decreasing)
    return Cp * (rec_rate - ion_rate) / ((1 + z) * Hz)


# ══════════════════════════════════════════════════════════════════════════
#  Solve
# ══════════════════════════════════════════════════════════════════════════

print("=" * 72)
print("  GERT Paper 3 — Ξ(α) with Peebles Recombination")
print("=" * 72)

z_start = 1800
z_end   = 50
Xe_init = float(Xe_saha(z_to_alpha(z_start)))

print(f"\n  Integrating Peebles ODE from z = {z_start} to z = {z_end}")
print(f"  Initial X_e (Saha at z={z_start}) = {Xe_init:.6f}")

z_grid   = np.linspace(z_start, z_end, 200000)
sol      = odeint(dXe_dz_peebles, Xe_init, z_grid,
                  rtol=1e-11, atol=1e-14, mxstep=100000)
Xe_peebles_arr = np.clip(sol[:, 0], 1e-30, 1.0)

print(f"  Done. Freeze-out X_e(z={z_end}) = {Xe_peebles_arr[-1]:.4e}")


def Xe_peebles(alpha):
    """Interpolated Peebles ionization fraction X_e(α)."""
    z = alpha_to_z(alpha)
    if z > z_start:
        return float(Xe_saha(alpha))
    if z < z_end:
        return float(Xe_peebles_arr[-1])
    return float(np.interp(z, z_grid[::-1], Xe_peebles_arr[::-1]))


# ══════════════════════════════════════════════════════════════════════════
#  Comparison table
# ══════════════════════════════════════════════════════════════════════════

print(f"\n{'='*72}")
print(f"  Ionization history: Saha vs Peebles")
print(f"{'='*72}")
print(f"\n  {'α':<10} {'z':<8} {'T(K)':<8} {'X_e(Saha)':<13} {'X_e(Peebles)':<14} {'P/S'}")
print(f"  {'-'*62}")

for z in [1800, 1500, 1400, 1300, 1200, 1150, 1100, 1080, 1050,
          1000, 900, 800, 600, 400, 200]:
    alpha = z_to_alpha(z)
    T  = T0 * (1 + z)
    xs = float(Xe_saha(alpha))
    xp = Xe_peebles(alpha)
    ps = xp / xs if xs > 1e-30 else float('inf')
    mark = " ← slower" if 1.5 < ps < 1000 else (" ← freeze-out" if ps > 1000 else "")
    print(f"  {alpha:<10.4f} {z:<8d} {T:<8.0f} {xs:<13.4e} {xp:<14.4e} {ps:<.2f}{mark}")


# ══════════════════════════════════════════════════════════════════════════
#  Ξ(α) profile
# ══════════════════════════════════════════════════════════════════════════

print(f"\n{'='*72}")
print(f"  Ξ(α) comparison")
print(f"{'='*72}")
print(f"\n  {'α':<10} {'z':<8} {'Ξ(Saha)':<14} {'Ξ(Peebles)':<14}")
print(f"  {'-'*48}")

for alpha in [-3.3, -3.2, -3.15, -3.10, -3.08, -3.06, -3.05, -3.04,
              -3.03, -3.02, -3.00, -2.95, -2.90, -2.85, -2.80, -2.50]:
    z    = alpha_to_z(alpha)
    xi_s = Xi(alpha, Xe_func=Xe_saha)
    xi_p = Xi(alpha, Xe_func=Xe_peebles)

    marks = []
    if 0.3 < xi_s < 3: marks.append("S")
    if 0.3 < xi_p < 3: marks.append("P")
    mark = f"  ← {','.join(marks)} ≈ 1" if marks else ""

    print(f"  {alpha:<10.4f} {z:<8.0f} {xi_s:<14.4e} {xi_p:<14.4e}{mark}")


# ══════════════════════════════════════════════════════════════════════════
#  Emergence boundaries
# ══════════════════════════════════════════════════════════════════════════

print(f"\n{'='*72}")
print(f"  Emergence boundary α_em")
print(f"{'='*72}")

results = {}
for label, Xe_f in [("Saha", Xe_saha), ("Peebles", Xe_peebles)]:
    try:
        a_em = brentq(lambda a: Xi(a, Xe_func=Xe_f) - 1.0,
                      -3.3, -2.5, xtol=1e-12)
        z_em = alpha_to_z(a_em)
        T_em = alpha_to_T(a_em)
        Xe_em = Xe_f(a_em)
        results[label] = dict(alpha=a_em, z=z_em, T=T_em, Xe=Xe_em)

        print(f"\n  {label}:")
        print(f"    α_em  = {a_em:.6f}")
        print(f"    z_em  = {z_em:.1f}")
        print(f"    T_em  = {T_em:.0f} K  ({T_em*k_B/eV:.4f} eV)")
        print(f"    X_e   = {Xe_em:.6f}")
    except Exception as e:
        print(f"\n  {label}: {e}")

# ── Combined result ─────────────────────────────────────────────────────
if len(results) == 2:
    a_s = results["Saha"]["alpha"]
    a_p = results["Peebles"]["alpha"]
    mid = (a_s + a_p) / 2
    half = abs(a_p - a_s) / 2

    print(f"\n{'='*72}")
    print(f"  Combined result")
    print(f"{'='*72}")
    print(f"\n  α_em (Saha)    = {a_s:.2f}")
    print(f"  α_em (Peebles) = {a_p:.2f}")
    print(f"  α_em (combined)= {mid:.1f} ± {half:.1f}")
    print(f"\n  Δα = {abs(a_p - a_s):.2f}  "
          f"({abs(a_p-a_s)/(alpha_crit-mid)*100:.1f}% of relativistic domain)")
    print(f"\n  Relativistic domain: [{mid-half:.1f}, {alpha_crit:.2f}]"
          f" = {alpha_crit-mid:.1f} ± {half:.1f} decades")
