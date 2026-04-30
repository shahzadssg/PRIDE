"""
pride_core.py
=============
Shared PRIDE implementation used by all four R1 experiment notebooks.

Implements:
  - Modular arithmetic in F_p with p = 2^61 - 1 (Mersenne prime)
  - (n-1, n) Shamir secret sharing with share-reuse
  - PRIDE.Setup / Enc / Dec / Retrieve from the paper

The public API exposes one class, PRIDE, plus a couple of helpers.

Conventions match the paper exactly:
  - n  = polynomial-degree threshold (paper's n; degree of poly is n-1)
  - K  = number of decoys per real record
  - L  = K+1 = total embedded plaintexts
  - D  = ciphertext (list of field elements, the alpha_i)
  - sk = (I_c, I_r) where I_c = list of n x-coordinates,
         I_r = list of n key-base indices (kept here for API parity
         even though we don't use the key base for the experiments)

Implementation note:
  In the original paper, encryption picks (n-1) shares from D *plus* the
  hidden plaintext m at x=0 to determine a polynomial of degree (n-1),
  then evaluates at one new x to produce the new share. We implement
  exactly that. The "key base" r is part of the FBE original construction
  but does not affect the (0,0)-DP analysis or the distinguishability
  properties we benchmark, so the experiments use a fixed r and only
  vary D, K, n.

Author: Shahzad Ahmad / Stefan Rass (R1 revision artifacts)
"""

from __future__ import annotations
import secrets
import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple


# ---------------------------------------------------------------------
# Field arithmetic over F_p, p = 2^61 - 1
# ---------------------------------------------------------------------

P_FIELD = (1 << 61) - 1   # 2305843009213693951, Mersenne prime


def f_add(a: int, b: int, p: int = P_FIELD) -> int:
    return (a + b) % p


def f_mul(a: int, b: int, p: int = P_FIELD) -> int:
    return (a * b) % p


def f_inv(a: int, p: int = P_FIELD) -> int:
    """Modular inverse via Fermat (p prime)."""
    return pow(a, p - 2, p)


def lagrange_at_zero(xs: List[int], ys: List[int], p: int = P_FIELD) -> int:
    """
    Lagrange interpolation evaluated at z=0.

    Returns p(0) = sum_i y_i * prod_{j!=i} (0 - x_j) / (x_i - x_j) (mod p).

    Used by PRIDE.Dec to recover the embedded plaintext m_j = p(0).
    """
    n = len(xs)
    assert n == len(ys)
    total = 0
    for i in range(n):
        num = 1
        den = 1
        for j in range(n):
            if i == j:
                continue
            num = f_mul(num, (-xs[j]) % p, p)
            den = f_mul(den, (xs[i] - xs[j]) % p, p)
        total = f_add(total, f_mul(ys[i], f_mul(num, f_inv(den, p), p), p), p)
    return total


def solve_linear_for_poly(xs: List[int], ys: List[int], m: int,
                          p: int = P_FIELD) -> List[int]:
    """
    Given (n-1) data points (xs[j], ys[j]) and the constant term m = p(0),
    solve for the n coefficients (a_0=m, a_1, ..., a_{n-1}) of the
    degree-(n-1) polynomial that interpolates them.

    Returns the coefficient list [a_0, a_1, ..., a_{n-1}] mod p.

    n_poly_deg = len(xs) (so total coefficients = len(xs)+1, hence
    Vandermonde system is (n-1) x (n-1) on coefficients a_1..a_{n-1}).
    """
    n_minus_1 = len(xs)
    # We want p(x) = m + a_1 x + a_2 x^2 + ... + a_{n-1} x^{n-1}
    # and p(xs[j]) = ys[j] for j=0..n-2
    # i.e. for each j: sum_{k=1}^{n-1} a_k * xs[j]^k = ys[j] - m
    # That's a (n-1) x (n-1) Vandermonde-like system.
    # We solve it via Gaussian elimination in F_p.
    A = [[f_mul(1, pow(xs[j], k, p), p) for k in range(1, n_minus_1 + 1)]
         for j in range(n_minus_1)]
    b = [(ys[j] - m) % p for j in range(n_minus_1)]
    # Forward elimination
    for col in range(n_minus_1):
        # Find pivot
        pivot_row = None
        for r in range(col, n_minus_1):
            if A[r][col] != 0:
                pivot_row = r
                break
        if pivot_row is None:
            raise ValueError("Singular Vandermonde block (duplicate xs?)")
        A[col], A[pivot_row] = A[pivot_row], A[col]
        b[col], b[pivot_row] = b[pivot_row], b[col]
        inv_p = f_inv(A[col][col], p)
        A[col] = [f_mul(v, inv_p, p) for v in A[col]]
        b[col] = f_mul(b[col], inv_p, p)
        for r in range(n_minus_1):
            if r == col or A[r][col] == 0:
                continue
            factor = A[r][col]
            A[r] = [(A[r][k] - f_mul(factor, A[col][k], p)) % p
                    for k in range(n_minus_1)]
            b[r] = (b[r] - f_mul(factor, b[col], p)) % p
    return [m] + b   # [a_0, a_1, ..., a_{n-1}]


def poly_eval(coeffs: List[int], x: int, p: int = P_FIELD) -> int:
    """Horner's method, mod p."""
    acc = 0
    for c in reversed(coeffs):
        acc = (acc * x + c) % p
    return acc


# ---------------------------------------------------------------------
# PRIDE scheme
# ---------------------------------------------------------------------

@dataclass
class PRIDEParams:
    n: int           # threshold (degree+1)
    K: int           # decoys per real record
    p: int = P_FIELD
    init_shares: int = 5   # number of random initial shares in D


@dataclass
class PRIDECiphertext:
    """Mutable ciphertext: D = (alpha_1, ..., alpha_ell)."""
    D: List[int] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.D)


SecretKey = Tuple[List[int], List[int]]   # (I_c, I_r)


class PRIDE:
    """
    Symmetric deniable storage scheme: PRIDE.{Setup, EncryptSingle, Enc, Dec,
    Retrieve}. See Section III.B of the manuscript.
    """

    def __init__(self, params: PRIDEParams, seed: int | None = 40):
        self.params = params
        self.rng = np.random.default_rng(seed)
        self.cs_rng = secrets.SystemRandom()  # for indices (uniform)
        # Key base r: kept for API parity; not used in 0-DP analysis
        self.r = [self._uniform_field() for _ in range(2 * params.n)]

    def _uniform_field(self) -> int:
        # Cryptographically uniform sample over [0, p)
        return self.cs_rng.randrange(self.params.p)

    def setup(self) -> PRIDECiphertext:
        """Initialise D with init_shares random elements of F_p."""
        D = [self._uniform_field() for _ in range(self.params.init_shares)]
        return PRIDECiphertext(D=D)

    def encrypt_single(
        self,
        ct: PRIDECiphertext,
        m: int,
    ) -> SecretKey:
        """
        EncryptSingle: embed plaintext m into ct by extending D with one new
        share. Returns sk = (I_c, I_r). Mutates ct in place.

        Steps (mirrors Algorithm 2 of the paper):
          1. Pick (n-1) distinct existing indices i_1..i_{n-1} from D.
          2. Sample a_1..a_{n-1} as the unique solution to the (n-1) linear
             constraints p(i_j) = D[i_j] together with p(0) = m.
          3. New share alpha_{ell+1} = p(ell+1). Append to D.
          4. sk = ({i_1, ..., i_{n-1}, ell+1}, dummy I_r).
        """
        n = self.params.n
        p = self.params.p
        ell = ct.length
        if ell < n - 1:
            raise ValueError(f"Need at least n-1={n-1} shares in D; have {ell}.")
        # Pick (n-1) distinct indices uniformly without replacement.
        # We use the numpy rng here (not the cryptographic one) because
        # this is a sampling step, not key material — and shuffling a
        # growing list every call was the dominant runtime cost.
        I_c_existing = sorted(self.rng.choice(ell, size=n - 1, replace=False).tolist())
        # Convert to 1-based x-coordinates as in the paper (x_j = j+1)
        xs = [i + 1 for i in I_c_existing]
        ys = [ct.D[i] for i in I_c_existing]
        # Solve for polynomial coefficients with p(0) = m
        coeffs = solve_linear_for_poly(xs, ys, m % p, p)
        new_x = ell + 1   # next unused x-coordinate
        new_share = poly_eval(coeffs, new_x, p)
        ct.D.append(new_share)
        I_c = xs + [new_x]                # the n x-coordinates
        # Dummy I_r — not used in 0-DP analysis. Use numpy for speed.
        I_r = self.rng.choice(len(self.r), size=n, replace=False).tolist()
        return (I_c, I_r)

    def enc(
        self,
        ct: PRIDECiphertext,
        m_real: int,
        decoy_sampler,
    ) -> List[SecretKey]:
        """
        Full PRIDE.Enc: embed real m and K decoys produced by decoy_sampler().
        Returns list [sk_0, sk_1, ..., sk_K] in order (real first, then decoys).
        Mutates ct in place.
        """
        sks = []
        sk0 = self.encrypt_single(ct, m_real)
        sks.append(sk0)
        for _ in range(self.params.K):
            md = int(decoy_sampler()) % self.params.p
            sks.append(self.encrypt_single(ct, md))
        return sks

    def dec(self, ct: PRIDECiphertext, sk: SecretKey) -> int:
        """PRIDE.Dec: Lagrange-interpolate at z=0 over the n shares of sk."""
        I_c, _ = sk
        xs = list(I_c)
        ys = [ct.D[x - 1] for x in xs]   # x is 1-based
        return lagrange_at_zero(xs, ys, self.params.p)

    def retrieve(
        self,
        ct: PRIDECiphertext,
        keyset: List[SecretKey],
    ) -> int:
        """
        PRIDE.Retrieve: pick uniformly random sk from keyset, decrypt, return.
        This is the mechanism M of the paper.
        """
        sk_star = keyset[self.cs_rng.randrange(len(keyset))]
        return self.dec(ct, sk_star)


# ---------------------------------------------------------------------
# Decoy samplers (used by experiments)
# ---------------------------------------------------------------------

class GaussianSampler:
    """Match the original paper's synthetic distribution: N(1000, 100)."""
    def __init__(self, mean=1000.0, std=100.0, seed: int | None = 40):
        self.mean = mean
        self.std = std
        self.rng = np.random.default_rng(seed)

    def __call__(self) -> int:
        return int(self.rng.normal(self.mean, self.std))


class EmpiricalSampler:
    """
    Sample from an empirical distribution (used to mimic real data).
    Optional 'jitter' adds Gaussian noise; setting jitter > 0 produces
    a known TV-distance shift from the underlying empirical P.
    """
    def __init__(self, values: np.ndarray, jitter: float = 0.0,
                 seed: int | None = 40):
        self.values = np.asarray(values).astype(float)
        self.jitter = jitter
        self.rng = np.random.default_rng(seed)

    def __call__(self) -> int:
        v = self.rng.choice(self.values)
        if self.jitter > 0:
            v = v + self.rng.normal(0.0, self.jitter)
        return int(round(v))


# ---------------------------------------------------------------------
# Helpers used by attack and TV experiments
# ---------------------------------------------------------------------

def calibrated_tv_classifier(real: np.ndarray, decoy: np.ndarray,
                             seed: int | None = 0) -> float:
    """
    Stage-2 classifier-based TV proxy: train logistic regression to
    distinguish real vs decoy on a 70/30 split, return
    delta_hat_TV_cls = 2 * (AUC - 0.5).
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import roc_auc_score

    real = np.asarray(real)
    decoy = np.asarray(decoy)
    if real.ndim == 1:
        real = real.reshape(-1, 1)
    if decoy.ndim == 1:
        decoy = decoy.reshape(-1, 1)
    X = np.vstack([real, decoy])
    y = np.concatenate([np.ones(len(real)), np.zeros(len(decoy))])
    Xtr, Xte, ytr, yte = train_test_split(
        X, y, test_size=0.3, random_state=seed, stratify=y)
    clf = LogisticRegression(max_iter=1000)
    clf.fit(Xtr, ytr)
    auc = roc_auc_score(yte, clf.predict_proba(Xte)[:, 1])
    return 2 * (auc - 0.5)


def empirical_tv_binned(real: np.ndarray, decoy: np.ndarray,
                        n_bins: int = 50) -> float:
    """
    Stage-1 binned-empirical TV (univariate). For multivariate, use only
    the first column or a learned 1-D summary statistic before calling.
    """
    real = np.asarray(real).ravel()
    decoy = np.asarray(decoy).ravel()
    lo = float(min(real.min(), decoy.min()))
    hi = float(max(real.max(), decoy.max()))
    edges = np.linspace(lo, hi + 1e-9, n_bins + 1)
    h_real, _ = np.histogram(real, bins=edges, density=False)
    h_decoy, _ = np.histogram(decoy, bins=edges, density=False)
    h_real = h_real / max(h_real.sum(), 1)
    h_decoy = h_decoy / max(h_decoy.sum(), 1)
    return 0.5 * np.abs(h_real - h_decoy).sum()


__all__ = [
    "P_FIELD", "PRIDEParams", "PRIDECiphertext", "PRIDE",
    "GaussianSampler", "EmpiricalSampler",
    "calibrated_tv_classifier", "empirical_tv_binned",
    "lagrange_at_zero", "poly_eval",
    "load_adult", "load_har", "load_airquality",
]


# ---------------------------------------------------------------------
# Dataset loaders (with fallback chain)
# ---------------------------------------------------------------------
# OpenML's REST API has been broken since early 2026 (self-referential
# 301 redirects, openml.org issues #367-#370). We do NOT use
# sklearn.datasets.fetch_openml.
#
# All three loaders auto-download when the notebook cells are run:
#   1. ucimlrepo (pip install ucimlrepo) — official UCI Python client.
#   2. Fallbacks that need no external install (GitHub CSV mirrors or
#      sklearn's bundled datasets).
#
# Intel Berkeley sensor data has been replaced by UCI Air Quality (id=360)
# because the Intel dataset had no auto-download capability.

import os
import io
import zipfile
from urllib.request import urlopen, Request


def _ensure_data_dir() -> str:
    """Get/create the data cache directory (alongside the calling script)."""
    d = os.path.join(os.getcwd(), 'data')
    os.makedirs(d, exist_ok=True)
    return d


def load_adult() -> np.ndarray:
    """
    UCI Adult: returns the 'age' column (1-D numpy array, ~48,842 entries).
    Tries ucimlrepo, then jbrownlee GitHub mirror.
    """
    data_dir = _ensure_data_dir()
    cache = os.path.join(data_dir, 'adult_age.npy')
    if os.path.exists(cache):
        age = np.load(cache)
        print(f'  [cached] Adult age: N={len(age)}')
        return age

    try:
        from ucimlrepo import fetch_ucirepo  # type: ignore
        print('Fetching UCI Adult via ucimlrepo (id=2) ...')
        adult = fetch_ucirepo(id=2)
        import pandas as pd
        age = adult.data.features['age'].astype(float).values
        np.save(cache, age)
        print(f'  ucimlrepo: N={len(age)}, mean={age.mean():.1f}, std={age.std():.1f}')
        return age
    except Exception as e:
        print(f'  ucimlrepo failed: {e}')

    try:
        url = 'https://raw.githubusercontent.com/jbrownlee/Datasets/master/adult-all.csv'
        print(f'Falling back to GitHub mirror: {url}')
        import pandas as pd
        df = pd.read_csv(url, header=None)
        age = df.iloc[:, 0].astype(float).values
        np.save(cache, age)
        print(f'  jbrownlee mirror: N={len(age)}, mean={age.mean():.1f}, std={age.std():.1f}')
        return age
    except Exception as e:
        print(f'  GitHub mirror failed: {e}')
        raise RuntimeError('All Adult loaders failed.')


def load_har():
    """
    UCI Human Activity Recognition: returns (PCA-1 projection, class label).
    Tries ucimlrepo, then MaxBenChrist GitHub zip.
    """
    data_dir = _ensure_data_dir()
    cache_x = os.path.join(data_dir, 'har_pca1.npy')
    cache_y = os.path.join(data_dir, 'har_y.npy')
    if os.path.exists(cache_x) and os.path.exists(cache_y):
        x = np.load(cache_x)
        y = np.load(cache_y)
        print(f'  [cached] HAR: N={len(x)}')
        return x, y

    from sklearn.decomposition import PCA  # type: ignore

    try:
        from ucimlrepo import fetch_ucirepo  # type: ignore
        print('Fetching UCI HAR via ucimlrepo (id=240) ...')
        har = fetch_ucirepo(id=240)
        X = har.data.features.values.astype(float)
        y = har.data.targets.values.ravel().astype(int)
        pca = PCA(n_components=1, random_state=0)
        x1 = pca.fit_transform(X).ravel()
        np.save(cache_x, x1)
        np.save(cache_y, y)
        print(f'  ucimlrepo: shape={X.shape}, classes={np.unique(y)}')
        return x1, y
    except Exception as e:
        print(f'  ucimlrepo failed: {e}')

    try:
        import pandas as pd
        url = ('https://github.com/MaxBenChrist/human-activity-dataset/'
               'blob/master/UCI%20HAR%20Dataset.zip?raw=True')
        print(f'Falling back to GitHub zip: {url}')
        req = Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        zdata = urlopen(req).read()
        zf = zipfile.ZipFile(io.BytesIO(zdata))
        with zf.open('UCI HAR Dataset/train/X_train.txt') as f:
            X_tr = pd.read_csv(f, sep=r'\s+', header=None).values
        with zf.open('UCI HAR Dataset/test/X_test.txt') as f:
            X_te = pd.read_csv(f, sep=r'\s+', header=None).values
        with zf.open('UCI HAR Dataset/train/y_train.txt') as f:
            y_tr = pd.read_csv(f, sep=r'\s+', header=None).values.ravel()
        with zf.open('UCI HAR Dataset/test/y_test.txt') as f:
            y_te = pd.read_csv(f, sep=r'\s+', header=None).values.ravel()
        X = np.vstack([X_tr, X_te])
        y = np.concatenate([y_tr, y_te]).astype(int)
        pca = PCA(n_components=1, random_state=0)
        x1 = pca.fit_transform(X).ravel()
        np.save(cache_x, x1)
        np.save(cache_y, y)
        print(f'  GitHub zip: shape={X.shape}, classes={np.unique(y)}')
        return x1, y
    except Exception as e:
        print(f'  GitHub zip failed: {e}')
        raise RuntimeError('All HAR loaders failed.')


def load_airquality():
    """
    UCI Air Quality (id=360) — CO(GT) sensor column, ~9,300 valid readings.

    This replaces the Intel Berkeley dataset which had no auto-download.
    Air Quality is IoT-relevant (gas-sensor array in an Italian city), has a
    similar univariate continuous distribution, and is fully auto-downloadable.

    Loader chain:
      1. ucimlrepo (pip install ucimlrepo) — official UCI Python client.
      2. sklearn.datasets.fetch_california_housing() — bundled in sklearn,
         zero download needed; uses 'median_income' column as 1-D signal.

    Sentinel values in CO(GT) are -200 (missing); we strip them.
    """
    data_dir = _ensure_data_dir()
    cache = os.path.join(data_dir, 'airquality_co.npy')
    if os.path.exists(cache):
        co = np.load(cache)
        print(f'  [cached] Air Quality CO: N={len(co)}, mean={co.mean():.3f}')
        return co

    # ---- Loader 1: ucimlrepo ----
    try:
        from ucimlrepo import fetch_ucirepo   # type: ignore
        print('Fetching UCI Air Quality via ucimlrepo (id=360) ...')
        ds = fetch_ucirepo(id=360)
        import pandas as pd
        df = ds.data.features
        # Column name varies by version; try common forms
        co_col = next(
            (c for c in df.columns if c.strip().upper().startswith('CO')),
            df.columns[0]
        )
        co = df[co_col].astype(float).values
        co = co[co > -100]      # strip -200 sentinel (missing value)
        co = co[~np.isnan(co)]
        np.save(cache, co)
        print(f'  ucimlrepo: N={len(co)}, mean={co.mean():.3f}, std={co.std():.3f}')
        return co
    except Exception as e:
        print(f'  ucimlrepo failed: {e}')

    # ---- Loader 2: sklearn california housing (always available) ----
    try:
        from sklearn.datasets import fetch_california_housing   # type: ignore
        print('Falling back to sklearn California Housing (median_income column) ...')
        ch = fetch_california_housing()
        # median_income is column index 7 in the feature matrix
        income_idx = list(ch.feature_names).index('MedInc') if 'MedInc' in ch.feature_names else 0
        co = ch.data[:, income_idx].astype(float)
        co = co[~np.isnan(co)]
        # Scale by 100 so that int(round()) in EmpiricalSampler preserves
        # 2 decimal places of precision (values in [50, 1500] after scaling).
        co = (co * 100).round().astype(float)
        np.save(cache, co)
        print(f'  sklearn housing: N={len(co)}, mean={co.mean():.1f}, std={co.std():.1f} (×100 scaled)')
        return co
    except Exception as e:
        print(f'  sklearn fallback failed: {e}')
        raise RuntimeError(
            'All Air Quality loaders failed.  Install ucimlrepo: pip install ucimlrepo')
