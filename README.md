# PRIDE: PRIvacy through DEniability

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![IEEE IoT-J](https://img.shields.io/badge/venue-IEEE%20Internet%20of%20Things%20Journal-blue)](https://ieee-iotj.org/)

Reference implementation for the paper:

> **PRIDE: PRIvacy through DEniability**  
> Shahzad Ahmad, Stefan Rass  
> LIT Secure and Correct Systems Lab, Johannes Kepler University Linz  
> *IEEE Internet of Things Journal* — Manuscript IoT-61935-2026 (under review)

PRIDE is a symmetric deniable storage scheme that combines **False-Bottom Encryption (FBE)** with **(n−1, n) polynomial secret sharing** over the Mersenne prime field $\mathbb{F}_{2^{61}-1}$. Any embedded plaintext can be opened with a different key, and the resulting view is statistically identical to opening any other plaintext — providing simulation-based deniability with a $(0,0)$-differential privacy corollary.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Core API](#core-api)
- [Experiment Notebooks](#experiment-notebooks)
- [Dataset Loaders](#dataset-loaders)
- [Reproducing the Paper's Results](#reproducing-the-papers-results)
- [Design Notes](#design-notes)
- [Citation](#citation)

---

## Overview

PRIDE stores $K+1$ plaintexts (one real, $K$ decoys) inside a single ciphertext $D$ — a list of field elements. Encrypting a new plaintext $m$ picks $n-1$ existing shares from $D$, solves for the unique degree-$(n-1)$ polynomial passing through those shares with $p(0) = m$, then appends one new evaluation point. Decryption with any of the $K+1$ secret keys recovers the corresponding plaintext via Lagrange interpolation at zero.

Key properties:

| Property | Value |
|---|---|
| Field | $\mathbb{F}_p$, $p = 2^{61}-1$ (Mersenne prime) |
| Secret sharing | $(n-1, n)$ Shamir with share reuse |
| Storage overhead | $K+1$ shares per real record |
| Encryption complexity | $O(n^3)$ per record (Gaussian elimination) |
| Retrieval complexity | $O(n^2)$ per key (Lagrange interpolation) |
| DP guarantee | $(\varepsilon=0, \delta=0)$-DP (corollary of simulation-based deniability) |

---

## Repository Structure

```
pride/
├── pride_core.py               # Core implementation (field arithmetic, PRIDE scheme, loaders)
│
├── param_sweep.ipynb        # Table IV  — parameter sweep over (n, K)
├── realdata.ipynb           # Table V   — calibrated TV proxy on real datasets
├── attack.ipynb             # Table VI  — decoy distinguishability attack
├── baselines.ipynb          # Table VII — comparison vs. vanilla FBE and Laplace-DP
├── performance_figures.ipynb# Figures   — scalability, encryption, decoy, storage plots
│
├── data/                       # Auto-created on first run; caches downloaded datasets
└── README.md
```

All notebooks are self-contained: datasets download automatically on first execution and are cached under `data/` for subsequent runs.

---

## Installation

Python 3.9 or later is required. Install dependencies with:

```bash
pip install numpy scikit-learn ucimlrepo pandas matplotlib
```

`ucimlrepo` is the official UCI ML Repository Python client. All dataset loaders have network-free fallbacks (GitHub mirrors or sklearn's bundled datasets), so the notebooks run even without a ucimlrepo connection.

---

## Quick Start

```python
from pride_core import PRIDE, PRIDEParams, GaussianSampler

# 1. Configure: threshold n=3, K=2 decoys per real record
params = PRIDEParams(n=3, K=2)
system = PRIDE(params, seed=40)
sampler = GaussianSampler(mean=1000.0, std=100.0, seed=40)

# 2. Initialise an empty ciphertext
ct = system.setup()

# 3. Encrypt a real message (m=42) alongside K=2 decoys
m_real = 42
secret_keys = system.enc(ct, m_real, sampler)
# secret_keys[0]  → key for the real plaintext
# secret_keys[1:] → keys for decoys (disclosed under coercion)

# 4. Decrypt with the real key → recovers m_real
recovered = system.dec(ct, secret_keys[0])
assert recovered == m_real % params.p

# 5. Decrypt with a decoy key → recovers a statistically indistinguishable value
decoy_val = system.dec(ct, secret_keys[1])
print(f"Real: {recovered},  Decoy: {decoy_val}")
```

---

## Core API

### `PRIDEParams`

```python
PRIDEParams(n, K, p=P_FIELD, init_shares=5)
```

| Parameter | Type | Description |
|---|---|---|
| `n` | `int` | Shamir threshold (polynomial degree + 1) |
| `K` | `int` | Number of decoys per real record |
| `p` | `int` | Field prime (default: $2^{61}-1$) |
| `init_shares` | `int` | Number of random shares pre-loaded into $D$ on setup |

### `PRIDE`

```python
system = PRIDE(params, seed=40)
```

| Method | Description |
|---|---|
| `setup()` | Initialise and return an empty `PRIDECiphertext` $D$ |
| `enc(ct, m, sampler)` | Embed real message `m` + `K` decoys into `ct`; return list of `K+1` secret keys |
| `encrypt_single(ct, m)` | Embed a single plaintext; return its secret key |
| `dec(ct, sk)` | Decrypt with secret key `sk`; return plaintext |
| `retrieve(ct, keys)` | Return all `K+1` plaintexts for a full key set |

### `GaussianSampler` / `EmpiricalSampler`

Decoy generators passed to `enc()`.

```python
# Exact distribution (for synthetic benchmarks)
sampler = GaussianSampler(mean=1000.0, std=100.0, seed=40)

# Empirical distribution from a dataset (with optional KDE jitter)
sampler = EmpiricalSampler(values=adult_age_array, jitter=0.0, seed=40)
```

### Field-level helpers

```python
lagrange_at_zero(xs, ys, p)     # Lagrange interpolation at z=0
poly_eval(coeffs, x, p)         # Horner's method evaluation
calibrated_tv_classifier(real, decoy, seed)   # δ̂_TV^cls = 2×(AUC − 0.5)
empirical_tv_binned(real, decoy, n_bins)      # Stage-1 binned TV proxy
```

---

## Experiment Notebooks

Each notebook is independent and produces one table from the paper. Run them in any order.

### `param_sweep.ipynb` → Table IV

Sweeps $n \in \{2,3,5,7,10\}$ and $K \in \{1,2,5,10,20\}$ at $N=10{,}000$ synthetic Gaussian records. Reports encryption time per record, retrieval time, and shares per record.

**Sanity check:** `shares_per_rec == K+1` for every cell.  
**Runtime:** ~5–10 minutes.

### `realdata.ipynb` → Table V

Encrypts $N=2{,}000$ records from each of three real datasets and measures the classifier-based TV proxy $\hat{\delta}_{\mathrm{TV}}^{\mathrm{cls}} = 2 \times (\mathrm{AUC} - 0.5)$.

**Sanity check:** Synthetic Gaussian AUC $\approx 0.500$.  
**Runtime:** ~2–4 minutes.

### `attack.ipynb` → Table VI

Runs 1,000 adversarial trials per (dataset, K) cell. The adversary trains a logistic-regression + random-forest ensemble on a held-out sample and tries to identify the real plaintext among $K+1$ candidates.

**Sanity check:** Synthetic Gaussian excess $\approx 0$ for all $K$.  
**Runtime:** ~10–20 minutes.

### `baselines.ipynb` → Table VII

Benchmarks PRIDE against vanilla FBE (no formal 0-DP claim) and Laplace-DP storage ($\varepsilon=0.1$, no deniability) on the UCI Adult dataset at $N=10{,}000$.

**Sanity check:** PRIDE / Laplace storage ratio $= K+1 = 3$.  
**Runtime:** ~2 minutes.

### `performance_figures.ipynb` → Figures

Generates the combined performance figure (scalability, encryption time, decoy generation time, storage) across $N \in \{1{,}000, \ldots, 100{,}000\}$ at $(n=3, K=2)$. Saves `pride_performance_combined.png`.

**Runtime:** ~10–15 minutes.

---

## Dataset Loaders

All three loaders are in `pride_core.py` and follow the same pattern: check the local cache (`data/`), try `ucimlrepo`, fall back to a network-free alternative.

| Function | Dataset | Primary | Fallback |
|---|---|---|---|
| `load_adult()` | UCI Adult — `age` column | `ucimlrepo` (id=2) | jbrownlee GitHub CSV |
| `load_har()` | UCI HAR — PCA-1 projection | `ucimlrepo` (id=240) | MaxBenChrist GitHub zip |
| `load_airquality()` | UCI Air Quality — CO sensor | `ucimlrepo` (id=360) | `sklearn` California Housing |

> **Note:** `fetch_openml` is not used anywhere. OpenML's REST API has had self-referential 301 redirect issues since early 2026.

Downloaded files are cached as `.npy` arrays under `data/` after the first run, so subsequent notebook executions are instant.

---

## Reproducing the Paper's Results

1. Clone the repository and install dependencies:
   ```bash
   git clone https://github.com/your-org/pride.git
   cd pride
   pip install numpy scikit-learn ucimlrepo pandas matplotlib
   ```

2. Open Jupyter and run the notebooks in any order:
   ```bash
   jupyter notebook
   ```

3. All four tables and the performance figure are fully reproducible from `pride_core.py` with `seed=40`. The `data/` directory is created automatically; no manual dataset placement is required.

**Reproducibility note:** Absolute timing numbers depend on hardware. Relative scaling behaviour (linear in $N$, cubic in $n$, linear in $K$) and all statistical results (AUC, $\hat{p}_{\mathrm{att}}$, TV proxy values) are seed-deterministic and hardware-independent.

---

## Design Notes

**Why $p = 2^{61}-1$?** This Mersenne prime supports fast modular reduction via bit-shift tricks and is compatible with Circle-STARK-based ZK proofs (which operate over the circle group of order $p+1 = 2^{61}$), making it a natural choice for the designated-verifier audit architecture described in the paper.

**Share reuse.** Each `enc()` call adds exactly one new share to $D$, regardless of $n$ or $K$. After encrypting $N$ real records with $K$ decoys each, `ct.length = init_shares + N*(K+1)`, giving exactly $K+1$ shares per record in the limit of large $N$.

**Non-cryptographic RNG for index sampling.** The `numpy` RNG is used only to select which $n-1$ existing shares to reuse. Actual share values and field elements use `secrets.SystemRandom`, which is cryptographically uniform over $[0, p)$.

**Why not `fetch_openml`?** OpenML's REST API has been unreliable since early 2026. The fallback chain (`ucimlrepo` → GitHub mirror → bundled sklearn data) ensures notebooks run without manual intervention.

---

## Citation

If you use this code in your research, please cite:

```bibtex
@article{ahmad_pride_2026,
  author  = {Ahmad, Shahzad and Rass, Stefan},
  title   = {{PRIDE}: {PRIvacy} through {DEniability}},
  journal = {IEEE Internet of Things Journal},
  year    = {2026},
  note    = {Manuscript IoT-61935-2026, under review}
}
```

---

## Authors

**Shahzad Ahmad** — LIT Secure and Correct Systems Lab, JKU Linz  
**Stefan Rass** — LIT Secure and Correct Systems Lab, JKU Linz

Questions about the implementation: open an issue or email `shahzad.ahmad@jku.at`.
