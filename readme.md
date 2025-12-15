
# Fractal-Cosmo

### A Self-Consistent Parametric Cosmology Framework Based on Fractal Holographic Duality

---

## Overview

**Fractal-Cosmo v17.1** is an open-source Python implementation of a  cosmological framework based on **Fractal Holographic Duality (FHD)**.  
The code accompanies the paper:

> *JWST Observational Anomalies and the Parametric Reconstruction of the Hubble Tension  
> Based on Fractal Holographic Duality *

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

### 4. Output

The code will generate:

- `output.txt`: Contains all numerical results, tables, and figures presented in the paper.

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


