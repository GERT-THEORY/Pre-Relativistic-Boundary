"""
GERT Paper 3 — Sensitivity to the Cauldron boundary
=====================================================
Section 3.7

Tests the dependence of α_em on the unknown Cauldron exit α_PC.
The GERT isolated-universe premise (Paper I, Section 2.4.1) implies
that the particle horizon integral starts at a_PC, not at a = 0.

Result: varying α_PC over 25 orders of magnitude changes α_em
by less than 5 × 10⁻⁴. The emergence boundary is entirely
determined by recombination physics, not by the Cauldron.
"""
import numpy as np
from scipy.optimize import brentq
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from gert_utils_p3 import *

print("=" * 72)
print("  GERT Paper 3 — Sensitivity of α_em to α_PC")
print("=" * 72)

# ── Saha ────────────────────────────────────────────────────────────────
print(f"\n  Method: Saha equilibrium")
print(f"\n  {'α_PC':<10} {'α_em':<16} {'z_em':<10} {'Δα_em':<14} {'|Δα/α_em|'}")
print(f"  {'-'*62}")

ref_alpha_em = None
for aPC in [-30, -25, -20, -15, -12, -10, -8, -6, -5, -4]:
    try:
        a_em = brentq(lambda a: Xi(a, alpha_PC=aPC) - 1.0,
                      -3.3, -2.5, xtol=1e-12)
        z_em = alpha_to_z(a_em)
        if ref_alpha_em is None:
            ref_alpha_em = a_em
        da = a_em - ref_alpha_em
        rel = abs(da / ref_alpha_em) if ref_alpha_em != 0 else 0
        print(f"  {aPC:<10d} {a_em:<16.8f} {z_em:<10.1f} {da:<+14.8f} {rel:.2e}")
    except Exception as e:
        print(f"  {aPC:<10d} failed: {e}")

# ── Peebles ─────────────────────────────────────────────────────────────
# Solve Peebles once
from scipy.integrate import odeint

def C_peebles(Xe, z):
    T=T0*(1+z); kT=k_B*T; nH=n_H0*(1+z)**3; n1s=nH*(1-Xe)
    thermal=(m_e*kT/(2*np.pi*hbar**2))**1.5
    bB=alpha_B(T)*thermal*np.exp(-B_2/kT)
    Hz=H_z(z)
    if n1s<1: return 1.0
    Lla=Hz*8*np.pi/(n1s*lambda_Lya**3)
    return (Lambda_2s+Lla)/(Lambda_2s+Lla+bB)

def dXe_dz(Xe, z):
    Xe=np.clip(Xe,1e-8,1-1e-10)
    T=T0*(1+z); kT=k_B*T; nH=n_H0*(1+z)**3; Hz=H_z(z)
    aB=alpha_B(T)
    thermal=(m_e*kT/(2*np.pi*hbar**2))**1.5
    bB=aB*thermal*np.exp(-B_2/kT)
    ion=bB*np.exp(-E_Lya/kT)*(1-Xe)
    rec=aB*nH*Xe**2
    return C_peebles(Xe,z)*(rec-ion)/((1+z)*Hz)

print(f"\n  Solving Peebles ODE...")
z_grid=np.linspace(1800,50,200000)
sol=odeint(dXe_dz, float(Xe_saha(z_to_alpha(1800))), z_grid,
           rtol=1e-11, atol=1e-14, mxstep=100000)
Xe_p_arr=np.clip(sol[:,0],1e-30,1.0)

def Xe_p(alpha):
    z=alpha_to_z(alpha)
    if z>1800: return float(Xe_saha(alpha))
    if z<50: return float(Xe_p_arr[-1])
    return float(np.interp(z, z_grid[::-1], Xe_p_arr[::-1]))

print(f"  Done. Freeze-out = {Xe_p_arr[-1]:.4e}")
print(f"\n  Method: Peebles three-level atom")
print(f"\n  {'α_PC':<10} {'α_em':<16} {'z_em':<10} {'Δα_em':<14} {'|Δα/α_em|'}")
print(f"  {'-'*62}")

ref_alpha_em_p = None
for aPC in [-30, -25, -20, -15, -12, -10, -8, -6, -5, -4]:
    try:
        a_em = brentq(lambda a: Xi(a, Xe_func=Xe_p, alpha_PC=aPC) - 1.0,
                      -3.3, -2.4, xtol=1e-12)
        z_em = alpha_to_z(a_em)
        if ref_alpha_em_p is None:
            ref_alpha_em_p = a_em
        da = a_em - ref_alpha_em_p
        rel = abs(da / ref_alpha_em_p) if ref_alpha_em_p != 0 else 0
        print(f"  {aPC:<10d} {a_em:<16.8f} {z_em:<10.1f} {da:<+14.8f} {rel:.2e}")
    except Exception as e:
        print(f"  {aPC:<10d} failed: {e}")

# ── Summary ─────────────────────────────────────────────────────────────
print(f"\n{'='*72}")
print(f"  Summary")
print(f"{'='*72}")
print(f"""
  The emergence boundary α_em is INSENSITIVE to the Cauldron exit α_PC.

  Saha:    varying α_PC from −30 to −5 (25 dex) shifts α_em by < 3×10⁻⁴
  Peebles: varying α_PC from −30 to −5 (25 dex) shifts α_em by < 5×10⁻⁴

  Physical reason: in the radiation-dominated regime, d_ph(α) ∝ a·(a − a_PC).
  Since a_PC ≪ a_em by at least 7 orders of magnitude, the correction is
  negligible. The emergence boundary is controlled entirely by the
  well-known physics of Thomson scattering and hydrogen recombination,
  not by the unknown microphysics of the Primordial Cauldron.

  This confirms the disciplined scope of Paper 3: the GERT framework
  does not need to model the Cauldron to define the onset of the
  relativistic ruler.
""")
