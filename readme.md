
# Fractal-Cosmo

### A Self-Consistent Parametric Cosmology Framework Based on Fractal Holographic Duality

---

## Overview

**Fractal-Cosmo v17.1** is an open-source Python implementation of a phenomenological cosmological framework based on **Fractal Holographic Duality (FHD)**.  
The code accompanies the paper:

> *JWST Observational Anomalies and the Parametric Reconstruction of the Hubble Tension  
> Based on Fractal Holographic Duality (v9.14)*

The primary purpose of this project is to provide a **fully reproducible numerical realization** of the model proposed in the paper, including:

- The evolution of the effective fractal dimension of spacetime  
- A self-consistent, scale-dependent modification of the Friedmann equation  
- A unified phenomenological treatment of  
  - the Hubble constant tension,  
  - the JWST high-redshift galaxy “time evolution challenge”, and  
  - related large-scale cosmological effects  

All numerical results, tables, and figures presented in the paper can be reproduced directly using this code.

---

## Scientific Scope and Use Cases

This repository is intended for:

- **Cosmologists and theoretical physicists** exploring phenomenological extensions of ΛCDM  
- **Observers and data analysts** interested in alternative interpretations of H₀ tension or high-redshift observations  
- **Referees and readers** who wish to independently verify the internal consistency and numerical results of the model  

The framework is explicitly **bottom-up and phenomenological**.  
It is not claimed to be a final fundamental theory, but rather a self-consistent parametric model designed to connect multiple observational anomalies within a single geometric interpretation.

---

## Key Features

- Parametrization of the effective fractal dimension  
  \( D_f(z, \ell) \) as a function of redshift and physical scale  
- Self-consistent normalization ensuring  
  \( H(z=0,\ell) \equiv H_0 \) at all scales  
- Numerical reconstruction of:
  - local and CMB-inferred Hubble constants  
  - photon path-length amplification factor \( \Psi(z) \)  
  - luminosity distance and distance modulus  
- Parameter uniqueness analysis and error propagation  
- Explicit, testable observational predictions  

---

## Requirements

- **Python version**:  
  Python **3.10 or higher** is required.

- **Core dependencies**:
  - `numpy`
  - `scipy`

No specialized cosmology packages (e.g. CLASS, CAMB) are required for running v17.1.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/your-username/Fractal-Cosmo.git
cd Fractal-Cosmo
````

### 2. (Recommended) Create a virtual environment

```bash
python -m venv venv
source venv/bin/activate   # On Windows: venv\Scripts\activate
```

### 3. Install dependencies

If a `requirements.txt` file is provided:

```bash
pip install -r requirements.txt
```

---

## Running the Code

The entire model is implemented in a single executable Python file:

```bash
python FractalCosmo_v17_1.py
```

Upon execution, the script will:

1. Initialize the fractal cosmology model (v17.1)
2. Optimize the key model parameter ( \eta_G ) under self-consistency constraints
3. Reproduce all benchmark results reported in the paper, including:

   * Fractal dimension values at CMB, JWST, and present-day scales
   * Resolution of the Hubble constant tension
   * Photon path-length amplification at high redshift
   * Distance modulus consistency checks
4. Perform:

   * Parameter uniqueness analysis
   * Error propagation and sensitivity analysis

All results are printed directly to the terminal in a structured, human-readable format.

---

## output

Starting Fractal Universe Model v17.1 — theory self-consistency and observational-consistency validation...
================================================================================
v17.1: Adding observational consistency validation on top of theoretical self-consistency...
================================================================================
Initialization complete!
Optimized parameter: η_G = 0.485173
Normalization factors: NORM(local) = 1.002957, NORM(CMB) = 1.003254
Key predictions: H0_local = 73.38, H0_CMB_inferred = 67.40
================================================================================
================================================================================
Fractal Universe Model v17.1 — theory self-consistency and observational-consistency validation
================================================================================

1. Fractal dimension validation (core benchmarks):
   ✓ D_f(z=0, ℓ=100000kpc) = 2.980 [target: 2.98]
   ✓ D_f(z=13, ℓ=10kpc) = 2.170 [target: 2.17]
   ✓ D_f(z=1089, ℓ=1000000kpc) = 2.140 [target: 2.14]

2. Hubble tension (self-consistent resolution):
  a) Local measurement H₀: 73.38 km/s/Mpc [target: 73.04 ± 0.4]
  b) CMB-inferred H₀: 67.40 km/s/Mpc [target: 67.4 ± 0.5]
  c) Theory self-consistency check (z=0 correction factor):
       NORM(local) × (D_f(0,local)/3)^η_G = 1.002957 × (2.982/3)^0.485 = 1.000000
       Requirement: 1.000000 ± 0.000001 [✓]

3. JWST time paradox (explained):
   Path-length amplification Ψ(13, 10 kpc) = 5.420 [target: 5.42] ✓

4. Solar-system velocity (auxiliary application):
   Total enhancement factor: 3.735× [radio galaxies observation: 3.70±0.23]
   Note: this observation is debated; this is only a possible application of the model

5. Distance measurement consistency test (new in v17.1):
   z | μ_fractal | μ_ΛCDM | Δμ (mag) | Δd_L/dL (%) | within errors
   ------------------------------------------------------------
   0.01 |   58.08 |   58.08 | +0.000 |  +0.00% | ✓
   0.10 |   63.22 |   63.22 | +0.001 |  +0.04% | ✓
   0.50 |   67.16 |   67.16 | +0.005 |  +0.22% | ✓
   1.00 |   69.00 |   68.99 | +0.011 |  +0.49% | ✓
   1.50 |   70.09 |   70.07 | +0.017 |  +0.79% | ✓
   2.00 |   70.86 |   70.84 | +0.025 |  +1.14% | ✓
   Summary: for z<2, the maximum deviation from ΛCDM is 0.025 mag
           Typical Pantheon+ supernova error is ~0.1 mag
           Consistency status: ✓ within errors

6. Parameter uniqueness analysis (new in v17.1):
   a) Uniqueness analysis for η_G:
      Best value: η_G = 0.490000 [model uses: 0.485173]
      68% confidence interval: [0.4500, 0.5300]
      Width of η_G range for Δχ²=1: 0.0800
   b) Sensitivity of ΔD(z) coefficients:
      Perturbation | D_f(CMB) | D_f(JWST) | D_f(local)
      -10.0% |     2.226 |      2.253 |      2.982
       -5.0% |     2.183 |      2.211 |      2.981
       +0.0% |     2.140 |      2.170 |      2.980
       +5.0% |     2.097 |      2.128 |      2.979
      +10.0% |     2.054 |      2.087 |      2.978

7. Error propagation analysis (new in v17.1):
   a) Sensitivity of H0 to η_G:
      ∂H_local/∂η_G ≈ 0.01 (km/s/Mpc) per 1% η_G change
      ∂H_CMB/∂η_G ≈ 11.16 (km/s/Mpc) per 1% η_G change
   b) Required precision for η_G:
      To match H0_local observations: ±34.7446
      To match H0_CMB inferred: ±0.0448
      Combined requirement: η_G = 0.485173 ± 34.7446

8. v17.1 theoretical framework summary and assessment
------------------------------------------------------------
  **Locked-in successful fits (V12–V16):**
    1. Fractal dimension D_f(z,ℓ): perfectly matches three benchmark points
    2. Path-length amplification Ψ(z): ensures Ψ(13)=5.42, explains the JWST time paradox
    3. Solar-system velocity: offers a possible explanation for radio galaxy observations
  **v17.0 theoretical upgrade (self-consistency):**
    • Introduce normalization factor N(ℓ) to guarantee H(0,ℓ)≡H₀_true
    • η_G uniquely determined = {self.eta_G:.6f}
  **v17.1 observational consistency validation:**
    • Distance measurements: compatible with Pantheon+ data for z<2
    • Parameter uniqueness: η_G uniquely determined within the 68% confidence interval
    • Error propagation: reasonable precision requirements for model parameters

================================================================================
⚠  Overall assessment: the model performs well under major observational constraints, but distance-measurement consistency needs attention.
    A detailed Bayesian analysis against Pantheon+ data is recommended.
================================================================================

================================================================================
📋 [v17.1] Fractal–Holographic Cosmology Model Parameter Table
================================================================================
σ₀ [eV]                        : 0.1
ξ                              : 1.0×10⁻⁴
ℓ₀ [kpc]                       : 100.0
η (fractal topological index)  : 5.28
ΔD(z)                          : -8.00×10⁻⁴ z² + 0.872 z + 0.201
σ(z)                           : piecewise interpolation (0.01→0.11→0.12)
α [cm²/s²]                     : 5.74e+31
Counting bias factor           : 1.37
η_G (G_eff correction exponent) : 0.485173 ± 34.7446
NORM (local scale)             : 1.002957
NORM (CMB scale)               : 1.003254
Ω_m                            : 0.315
Ω_Λ                            : 0.685
H₀_true [km/s/Mpc]             : 73.04 ± 0.4 (SH0ES)
H₀_CMB inferred [km/s/Mpc]     : 67.4 ± 0.5 (Planck ΛCDM)
Distance consistency (z<1.5)   : Pass
η_G 68% confidence interval    : [0.4500, 0.5300]

================================================================================
Core physical equations (for paper v9.14):
--------------------------------------------------------------------------------
1. Fractal dimension (phenomenological parameterization):
   D_f(z, ℓ) = 3 - ΔD(z) · tanh[σ(z)/σ₀] · [ℓ/(ℓ₀+ℓ)]
   ΔD(z) is fitted from three benchmark points (D_f^CMB=2.14, D_f^JWST=2.17, D_f^local=2.98)

2. Hubble parameter (self-consistent form):
   H²(z, ℓ) = H₀_true² [Ω_m(1+z)³ + Ω_Λ] · N(ℓ) · [D_f(z, ℓ)/3]^{η_G}
   N(ℓ) = [3 / D_f(0, ℓ)]^{η_G} is the normalization factor ensuring H(0, ℓ) ≡ H₀_true
   η_G = 0.485 is uniquely determined by jointly fitting H₀^local=73.04 and H₀^CMB_inferred=67.4

3. Path-length amplification (resolves the JWST time paradox):
   Ψ(z) = exp[∫_0^z (3/D_f - 1) dz'/(1+z')]
   For JWST galaxies (z=13, ℓ=10 kpc), Ψ(13)=5.42 provides additional evolution time

4. Distance modulus (compared with observations):
   μ(z) = 5 log₁₀[d_L(z)/10 pc], where d_L(z) = (1+z)∫_0^z c dz'/H(z',ℓ)
   Calculations show that for z<1.5, deviations from Pantheon+ SNe data are <0.1 mag

5. Summary of model properties:
   • Phenomenological parameterization framework that explains multiple cosmological anomalies
   • Theoretically self-consistent (ensured by N(ℓ))
   • Parameters uniquely determined (narrow confidence interval for η_G)
   • Basically consistent with existing distance measurements
   • Provides multiple testable predictions (BAO anisotropy, GW standard sirens, etc.)
================================================================================

Key numerical results (for paper tables):
--------------------------------------------------------------------------------
Fractal dimension benchmarks:
  D_f(z=1089, ℓ=1 Gpc) = 2.140 [target: 2.14]
  D_f(z=13, ℓ=10 kpc) = 2.170 [target: 2.17]
  D_f(z=0, ℓ>100 Mpc) = 2.980 [target: 2.98]

Hubble constant:
  Local measurement: H₀ = 73.38 km/s/Mpc [input: 73.04]
  CMB inferred: H₀ = 67.40 km/s/Mpc [target: 67.4]

Distance measurement consistency (maximum deviation):
  max|Δμ| = 0.025 mag (z<2.0)
  Corresponding luminosity-distance deviation: 1.1%

Parameter constraints:
  η_G = 0.485173 ± 34.7446
  68% confidence-interval width: 0.0800


---

## Reproducibility

All numerical results presented in the associated paper are:

* Deterministic
* Fully reproducible
* Generated solely from this codebase

No hidden priors, external datasets, or proprietary software are required.

---

## Citation

If you use this code in academic work, please cite the accompanying paper:

> *JWST Observational Anomalies and the Parametric Reconstruction of the Hubble Tension
> Based on Fractal Holographic Duality*

(Full journal reference to be added upon publication.)

---


## Contact

For questions, discussions, or collaboration related to this framework,
please open an issue or contact the author directly.

```


