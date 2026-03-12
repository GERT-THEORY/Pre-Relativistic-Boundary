# Pre-Relativistic-Boundary
Scripts relacionados ao paper  The Onset of the Relativistic Ruler: Metric Emergence and the Pre-Relativistic Boundary of the GERT Universe
# GERT Paper 3 — Supplementary Calculation Scripts

Supplementary numerical verification for:

**Dutra (2025). "The Onset of the Relativistic Ruler: Metric Emergence and the Pre-Relativistic Boundary of the GERT Universe."**

Paper III of the GERT series.

---

## Contents

| File                  | Sections | What it computes                                             |
| --------------------- | -------- | ------------------------------------------------------------ |
| `gert_utils_p3.py`    | —        | Shared constants, Paper 1 parameters, and common functions (Hubble parameter, Saha ionization, photon mean free path, particle horizon, emergence parameter $\Xi$) |
| `calc_Xi_saha.py`     | 3.5–3.6  | $\Xi(\alpha)$ profile with Saha equilibrium; emergence boundary $\alpha_{\rm em} = -3.047$ |
| `calc_Xi_peebles.py`  | 3.6      | Peebles three-level atom recombination; $\alpha_{\rm em} = -2.937$; combined result $\alpha_{\rm em} = -3.0 \pm 0.1$ |
| `calc_sensitivity.py` | 3.7      | Sensitivity of $\alpha_{\rm em}$ to the Cauldron boundary $\alpha_{\rm PC}$: $ |
| `calc_domain_map.py`  | 3.8–3.9  | Complete GERT domain map in natural time $\alpha = \log_{10}(a)$; structural symmetry with Paper II |

All scripts reproduce the numerical results of the paper exactly from the parameters of Paper 1.

---

## Requirements

- Python 3.7+
- NumPy
- SciPy

```bash
pip install numpy scipy
```

---

## Usage

All scripts import `gert_utils_p3.py` and must be run from the same directory:

```bash
python calc_Xi_saha.py
python calc_Xi_peebles.py
python calc_sensitivity.py
python calc_domain_map.py
```

---

## Key Parameters (from Paper 1)

| Parameter            | Value                      | Status                     |
| -------------------- | -------------------------- | -------------------------- |
| $H_0$                | 72.5 km/s/Mpc              | Empirical                  |
| $\Omega_{m,0}$       | 0.30                       | Fixed                      |
| $\Omega_{\Lambda,0}$ | 0.70                       | Fixed                      |
| $\Omega_{b,0}$       | 0.045                      | Fixed                      |
| $T_0$                | 2.725 K                    | Empirical                  |
| $\sigma_T$           | $6.652 \times 10^{-29}$ m² | QED exact                  |
| $B_H$                | 13.6 eV                    | Hydrogen ionization energy |

The uncertainty in $\alpha_{\rm em}$ is dominated by the treatment of recombination (Saha vs. Peebles).

---

## Key Results

### Emergence boundary

$$\alpha_{\rm em} = -3.0 \pm 0.1$$

| Method           | $\alpha_{\rm em}$ | $z_{\rm em}$ | $T_{\rm em}$ (K) | $X_e$ at emergence |
| ---------------- | ----------------- | ------------ | ---------------- | ------------------ |
| Saha equilibrium | $-3.047$          | 1114         | 3037             | 0.005              |
| Peebles 3-level  | $-2.937$          | 865          | 2358             | 0.008              |

### Insensitivity to the Primordial Cauldron

| $\alpha_{\rm PC}$ | $\Delta\alpha_{\rm em}$ (Saha) | $\Delta\alpha_{\rm em}$ (Peebles) |
| ----------------- | ------------------------------ | --------------------------------- |
| $-30$             | $0$                            | $0$                               |
| $-20$             | $0$                            | $0$                               |
| $-15$             | $0$                            | $0$                               |
| $-10$             | $< 10^{-6}$                    | $< 10^{-6}$                       |
| $-8$              | $< 10^{-6}$                    | $< 10^{-6}$                       |
| $-6$              | $3 \times 10^{-5}$             | $5 \times 10^{-5}$                |
| $-5$              | $3 \times 10^{-4}$             | $5 \times 10^{-4}$                |

### Complete GERT domain map

$$\underbrace{\alpha < \alpha_{\rm PC}}_{\text{Cauldron}} \quad \underbrace{\alpha_{\rm PC} < \alpha < \alpha_{\rm em}}_{\text{Pre-legible}} \quad \underbrace{\alpha_{\rm em} \leq \alpha \leq \alpha_{\rm crit}}_{\text{Relativistic regime}} \quad \underbrace{\alpha > \alpha_{\rm crit}}_{\text{Quasi-Vacuum}}$$

| Boundary            | $\alpha$ | Uncertainty   | Source                  |
| ------------------- | -------- | ------------- | ----------------------- |
| $\alpha_{\rm em}$   | $-3.0$   | $\pm 0.1$     | Recombination kinetics  |
| $\alpha_{\rm crit}$ | $+12.88$ | $\pm 0.12$    | $k_{\rm gas}$ (Paper I) |
| **Span**            | **15.9** | **$\pm$ 0.2** | **decades in $\alpha$** |

---

## Physical Definitions

### Metric emergence parameter

$$\Xi(\alpha) \equiv \frac{\lambda_\gamma(\alpha)}{d_{\rm ph}(\alpha)}$$

where $\lambda_\gamma$ is the photon mean free path and $d_{\rm ph}$ is the GERT particle horizon (with lower limit $\alpha_{\rm PC}$, not $-\infty$).

### GERT particle horizon

$$d_{\rm ph}(\alpha) = a \int_{a_{\rm PC}}^{a} \frac{c \, da'}{a'^2 \, H(\alpha')}$$

The lower limit reflects the isolated-universe premise (Paper I, Section 2.4.1): time is Work, and below $\alpha_{\rm PC}$ the integral has no ontological meaning. This is independently supported by the no-boundary proposal of Hartle & Hawking.

### Peebles C-factor

$$C = \frac{\Lambda_{2s} + \Lambda_{\rm Ly\alpha}}{\Lambda_{2s} + \Lambda_{\rm Ly\alpha} + \beta_B}$$

The C-factor captures the kinetic bottleneck: most atoms reaching $n=2$ are re-photoionized before decaying to the ground state ($C \sim 10^{-2}$ during recombination), slowing the process and producing a freeze-out residual $X_e \sim 3 \times 10^{-4}$.

---

## Technical Notes

### Peebles ODE formulation

The Peebles equation is integrated in redshift $z$ (decreasing):

$$\frac{dX_e}{dz} = \frac{C}{(1+z)H} \left[\alpha_B \, n_H \, X_e^2 - \beta_B \, e^{-E_{\rm Ly\alpha}/k_B T} \, (1 - X_e)\right]$$

The factor $e^{-E_{\rm Ly\alpha}/k_B T}$ ensures detailed balance with the Saha equation, since $B_2 + E_{\rm Ly\alpha} = B_1$. The solver uses `scipy.integrate.odeint` with `rtol=1e-11`, `atol=1e-14`.

### Relation to Paper 2 scripts

Paper 2 scripts use `gert_utils.py` with the gas-dominated Hubble parameter $H_{\rm gas}(a)$ for the ultra-dilute regime. Paper 3 scripts use `gert_utils_p3.py` with the standard radiation + matter + $\Lambda$ Hubble parameter, which is appropriate for the early Universe where the GERT $f_M$ and $f_L$ modifications are at their high-density plateaus.

---

## Related

- **Paper 1**: Gibbs Energy Redistribution Theory (GERT): A Thermodynamically Motivated Expansion History and the Hubble Tension
- **Paper 2**: Black Hole Thermodynamic Inversion and the Terminal State of the GERT Universe  
  DOI: [10.20944/preprints202603.0279.v1](https://doi.org/10.20944/preprints202603.0279.v1)
- **Paper 2 scripts**: [github.com/...] *(same repository, `/paper2/` directory)*
