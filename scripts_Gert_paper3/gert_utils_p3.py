"""
GERT Paper 3 — Shared utilities
=================================
Common constants, parameters, and functions used across all Paper 3 scripts.

Extends the Paper 2 utilities (gert_utils.py) with early-universe physics:
photon mean free path, particle horizon, ionization fraction, and the
metric emergence parameter Ξ(α).

All Paper 1 parameters are inherited. The GERT natural time coordinate
is α = log₁₀(a), with the thermodynamic control variable x = log₁₀(ρ).

IMPORTANT: This file uses the same Paper 1 best-fit parameters as
gert_utils.py from Paper 2. Any change to those parameters must be
propagated to both files.
"""

import numpy as np
from scipy.integrate import quad

# ══════════════════════════════════════════════════════════════════════════
#  Physical constants (SI)
# ══════════════════════════════════════════════════════════════════════════
G       = 6.674e-11        # gravitational constant [m³ kg⁻¹ s⁻²]
hbar    = 1.055e-34        # reduced Planck constant [J·s]
c       = 2.998e8          # speed of light [m/s]
k_B     = 1.381e-23        # Boltzmann constant [J/K]
m_e     = 9.109e-31        # electron mass [kg]
m_H     = 1.673e-27        # proton mass [kg]
M_sun   = 1.989e30         # solar mass [kg]
Mpc     = 3.085677581e22   # megaparsec [m]
yr      = 3.156e7          # year [s]
eV      = 1.602e-19        # electronvolt [J]
sigma_T = 6.652e-29        # Thomson cross section [m²]

# Hydrogen atomic physics
B_1       = 13.598 * eV    # ground-state ionization energy [J]
B_2       = B_1 / 4        # n=2 ionization energy [J]
E_Lya     = 3 * B_1 / 4    # Lyman-α transition energy [J]
Lambda_2s = 8.2245         # 2s→1s two-photon decay rate [s⁻¹]
lambda_Lya = 2 * np.pi * hbar * c / E_Lya  # Lyman-α wavelength [m]

# ══════════════════════════════════════════════════════════════════════════
#  Paper 1 best-fit parameters
# ══════════════════════════════════════════════════════════════════════════
H0_km      = 72.5          # Hubble constant [km/s/Mpc] — Empirical
Omega_L0   = 0.70          # dark energy density fraction — Fixed
Omega_m0   = 0.30          # matter density fraction — Fixed
Omega_b0   = 0.045         # baryon density fraction — Fixed
f_Lm       = 1.12          # entropic factor mid-value — Fixed
k_gas      = 0.143         # gas regime intensity — Free (+0.102/-0.103)
gamma_gas  = 0.50          # gas regime slope — Fixed
x_gas      = -26.750       # gas regime onset [log₁₀(ρ)] — Free (+0.219/-0.180)
T0         = 2.725         # CMB temperature today [K] — Empirical

# ══════════════════════════════════════════════════════════════════════════
#  Derived quantities
# ══════════════════════════════════════════════════════════════════════════
h           = H0_km / 100
H0          = H0_km * 1e3 / Mpc                         # [s⁻¹]
Omega_gamma = 2.469e-5 / h**2                            # photon density
Omega_nu    = Omega_gamma * 3.044 * (7/8) * (4/11)**(4/3)  # neutrino density
Omega_r0    = Omega_gamma + Omega_nu                     # total radiation
rho_crit    = 3 * H0**2 / (8 * np.pi * G)               # critical density [kg/m³]
rho_m0      = Omega_m0 * rho_crit                        # matter density today [kg/m³]
rho_b0      = Omega_b0 * rho_crit                        # baryon density today [kg/m³]
log_rho_m0  = np.log10(rho_m0)                           # log₁₀(ρ_m,0)

# Baryon/hydrogen number densities
n_b0        = rho_b0 / m_H                               # baryon number density [m⁻³]
n_H0        = n_b0 * 0.76                                # hydrogen number density [m⁻³]
                                                          # (He mass fraction Y ≈ 0.24)

# Characteristic photon scale
lambda_0    = 1.0 / (n_b0 * sigma_T)                     # [m]

# Paper 2 result
alpha_crit  = 12.8754                                     # ± 0.12
a_crit      = 10**alpha_crit


# ══════════════════════════════════════════════════════════════════════════
#  Coordinate conversions
# ══════════════════════════════════════════════════════════════════════════

def alpha_to_z(alpha):
    """Redshift z from GERT time α = log₁₀(a)."""
    return 10**(-alpha) - 1

def z_to_alpha(z):
    """GERT time α = log₁₀(a) from redshift z."""
    return -np.log10(1 + z)

def alpha_to_x(alpha):
    """GERT density variable x = log₁₀(ρ) from α."""
    return log_rho_m0 - 3 * alpha

def alpha_to_T(alpha):
    """CMB temperature [K] at GERT time α."""
    return T0 * 10**(-alpha)


# ══════════════════════════════════════════════════════════════════════════
#  Hubble parameter
# ══════════════════════════════════════════════════════════════════════════

def H_GERT(alpha):
    """GERT Hubble parameter H(α) [s⁻¹].

    For the early Universe (α ≪ 0, high density), the GERT fM and fL
    are at their high-density plateau values, and the standard
    radiation + matter + Λ form applies. The gas term and phase
    transitions in fM/fL activate only at lower densities (α > ~−1).
    """
    a = 10**alpha
    return H0 * np.sqrt(Omega_r0 / a**4 + Omega_m0 / a**3 + Omega_L0)

def H_z(z):
    """Hubble parameter H(z) [s⁻¹]."""
    return H0 * np.sqrt(Omega_r0*(1+z)**4 + Omega_m0*(1+z)**3 + Omega_L0)


# ══════════════════════════════════════════════════════════════════════════
#  Recombination physics
# ══════════════════════════════════════════════════════════════════════════

def alpha_B(T):
    """Case-B recombination coefficient α_B(T) [m³/s].

    Fitting formula from Péquignot, Petitjean & Boisson (1991),
    with Seager et al. (1999) fudge factor F = 1.14.
    Valid for 1000 K ≲ T ≲ 30000 K.
    """
    t4 = T / 1e4
    F = 1.14
    return F * 1e-19 * (4.309 / t4**0.6166) / (1 + 0.6703 * t4**0.5300)


def Xe_saha(alpha):
    """Saha equilibrium ionization fraction X_e(α).

    Solves X_e²/(1 − X_e) = S for hydrogen, where S is the Saha RHS.
    Returns X_e clipped to [1e−30, 1].
    """
    T  = alpha_to_T(alpha)
    kT = k_B * T
    nb = n_b0 * 10**(-3 * alpha)
    S  = (m_e * kT / (2 * np.pi * hbar**2))**1.5 * np.exp(-B_1 / kT) / nb
    Xe = (-S + np.sqrt(S**2 + 4*S)) / 2
    return np.clip(Xe, 1e-30, 1.0)


# ══════════════════════════════════════════════════════════════════════════
#  Photon mean free path and particle horizon
# ══════════════════════════════════════════════════════════════════════════

def lambda_gamma(alpha, Xe_func=None):
    """Photon mean free path λ_γ(α) [m].

    Parameters
    ----------
    alpha : float
        GERT time coordinate.
    Xe_func : callable, optional
        Ionization fraction function Xe(α). Defaults to Saha.
    """
    if Xe_func is None:
        Xe_func = Xe_saha
    Xe = max(Xe_func(alpha), 1e-30)
    a  = 10**alpha
    n_e = n_b0 / a**3 * Xe
    return 1.0 / (n_e * sigma_T)


def _conformal_integrand(a_prime):
    """Integrand for conformal time: 1/(a'² H(a'))."""
    alpha_prime = np.log10(a_prime)
    return 1.0 / (a_prime**2 * H_GERT(alpha_prime))


def d_ph(alpha, alpha_PC=-15):
    """GERT particle horizon d_ph(α) [m].

    The lower integration limit is a_PC = 10^{α_PC}, reflecting the
    GERT isolated-universe premise: the integral has no ontological
    meaning below the Cauldron exit.

    Parameters
    ----------
    alpha : float
        GERT time coordinate.
    alpha_PC : float
        Cauldron exit boundary in α (default: −15).
    """
    a    = 10**alpha
    a_PC = 10**alpha_PC
    if a <= a_PC:
        return 0.0
    eta, _ = quad(_conformal_integrand, a_PC, a, limit=300, epsrel=1e-10)
    return a * c * eta


# ══════════════════════════════════════════════════════════════════════════
#  Metric emergence parameter
# ══════════════════════════════════════════════════════════════════════════

def Xi(alpha, Xe_func=None, alpha_PC=-15):
    """Metric emergence parameter Ξ(α) = λ_γ(α) / d_ph(α).

    Ξ < 1 : pre-legible geometry (photons trapped)
    Ξ = 1 : emergence boundary α_em
    Ξ > 1 : relativistic effective regime (metric globally readable)

    Parameters
    ----------
    alpha : float
        GERT time coordinate.
    Xe_func : callable, optional
        Ionization fraction function Xe(α). Defaults to Saha.
    alpha_PC : float
        Cauldron boundary (default: −15).
    """
    lam = lambda_gamma(alpha, Xe_func)
    dph = d_ph(alpha, alpha_PC)
    if dph <= 0:
        return 0.0
    return lam / dph
