# PRIDE: PRIvacy through DEniability

## Overview

PRIDE is a framework that integrates **deniable encryption** with **0-differential privacy (0-DP)** to provide robust privacy guarantees for data stored in untrusted or coerced environments. This Python implementation demonstrates the core algorithms, including a modified False-Bottom Encryption (FBE) scheme based on polynomial secret sharing, and analyses their scalability and efficiency.

The goal is to enable data owners to achieve **plausible deniability**, where they can convincingly deny the true content of their stored data, even under duress or in the event of a breach. Simultaneously, it offers **perfect 0-DP**, ensuring that the data retrieval mechanism's output distribution is statistically identical regardless of the presence or absence of a specific user's real data.

### Key Features:

*   **0-Differential Privacy (0-DP):** Provides information-theoretic guarantees that retrieved data is perfectly private with respect to the true data, making it statistically indistinguishable from decoy data.
*   **Plausible Deniability & Coercion Resistance:** Allows users/providers to disclose "fake" plaintexts convincingly, protecting against physical or legal coercion.
*   **False-Bottom Encryption (FBE) Foundation:** Leverages a symmetric deniable encryption scheme to embed one real message and multiple decoy messages within a single ciphertext.
*   **Polynomial Secret Sharing:** Utilizes an \((n-1)\)-out-of-\(n\) Shamir-like secret sharing scheme for efficient share management and message reconstruction. The message `m` is embedded as `p(0)` of a polynomial of degree `n-1`.
*   **Post-Quantum Security Basis:** Built on information-theoretic secret sharing, which inherently provides resistance against quantum computing advancements.
*   **Benchmarked Scalability:** Experimental evaluation demonstrates linear scaling for encryption time, decoy generation time, and ciphertext storage with increasing dataset sizes, crucial for practical deployment.

## Problem Solved

Traditional encryption alone is insufficient against powerful adversaries who can compel users to reveal encryption keys. PRIDE addresses the urgent need for:

1.  **Coercion Resistance:** Empowering users to plausibly deny the true data even when forced to reveal decryption keys.
2.  **Leakage Resilience:** Protecting sensitive data from identification even if storage systems are compromised, by making the real data indistinguishable from decoys.
3.  **Accountability with Privacy:** Providing a framework that, while enabling deniability, also necessitates transparency and auditing mechanisms to prevent malicious data processors from abusing the 0-DP guarantee.

## Installation

This project is implemented in Python.

1.  **Clone the repository:**

    ```bash
    git clone https://github.com/your-username/PRIDE.git
    cd PRIDE
    ```

2.  **Install dependencies:**

    ```bash
    pip install numpy matplotlib
    ```

## Usage

The `PRIDE` class encapsulates the framework's functionalities. The main script (`pride_main.py` if you save it as such, or directly running the provided code) demonstrates setup, encryption, decoy generation, decryption, and benchmarks.

### Running the Example and Benchmarks

To run the full demonstration, including benchmarks and plot generation, execute the Python script:

```bash
python your_pride_script_name.py
```

Upon execution, you will see console output detailing the setup, benchmark progress, and retrieval demonstration. Performance plots (total processing time, encryption time, decoy generation time, and storage efficiency) will be saved in the `pride_plots/` directory.

### Code Structure

*   `mod_add(a, b)`, `mod_mul(a, b)`, `mod_inv(a)`: Modular arithmetic operations within the finite field \( \mathbb{Z}_p \).
*   `lagrange_interpolate(x, xs, ys)`: Implements Lagrange interpolation to reconstruct a polynomial and evaluate it at a given point `x`.
*   `solve_linear_system(matrix, vector, modulus)`: A helper function for Gaussian elimination to solve linear systems modulo `p`, used to find polynomial coefficients.
*   `PRIDE` Class:
    *   `__init__(self, n, K)`: Initializes the PRIDE instance with polynomial threshold `n` and decoy count `K`, and calls `setup()`.
    *   `setup(self)`: Initializes system parameters, the key base `r`, and the initial ciphertext `D`.
    *   `_encrypt_single(self, m)`: Internal method to encrypt a single message `m` using polynomial sharing, generating a secret key and updating the ciphertext `D`. This is where `p(0)=m` is established and coefficients for `p(x)` are determined.
    *   `PRIDE_enc(self, m, decoy_sampler)`: Encrypts a real message `m` and `K` decoy messages, collecting all generated secret keys and timing the process.
    *   `PRIDE_dec(self, secret_key)`: Decrypts a message using a provided secret key, reconstructing `p(0)` via Lagrange interpolation from `n` shares.
    *   `retrieval_mechanism(self)`: Simulates the data retrieval, randomly selecting a secret key (real or decoy) and decrypting it to demonstrate deniability.
*   `decoy_sampler_normal(mean, std)`: Samples a decoy message from a normal distribution.
*   `real_message_sampler(mean, std)`: Samples a real message from a similar distribution, consistent with the 0-DP assumption.

## Core Concepts Explained

### Deniable Encryption

Deniable encryption allows a party to plausibly deny the content of a plaintext that was encrypted, by being able to produce a different (fake) plaintext that also decrypts from the same ciphertext. PRIDE uses a *plan-ahead* deniable scheme, meaning decoys are prepared at the time of encryption.

### Differential Privacy (DP)

Differential Privacy offers a mathematical guarantee that queries on a dataset do not reveal specific information about any single individual. PRIDE achieves **0-Differential Privacy (0-DP)**, the strongest form of DP, ensuring that the act of data retrieval itself is perfectly private. This means the output distribution is identical whether or not a specific user's data is truly present.

### False-Bottom Encryption (FBE)

FBE is a symmetric deniable encryption technique where one ciphertext can simultaneously conceal multiple plaintexts (one real, several decoys). The core security of FBE relies on information-theoretic principles from secret sharing. PRIDE adapts FBE to integrate with polynomial secret sharing for its deniability and DP properties.

### Polynomial Secret Sharing (Shamir's Scheme)

PRIDE utilizes an \((n-1)\)-out-of-\(n\) Shamir secret sharing scheme. A secret is embedded as the constant term `p(0)` of a polynomial \(p(x)\) of degree \(n-1\). Any \(n\) points on this polynomial can reconstruct the secret, while \(n-1\) points reveal no information. This is fundamental to how multiple messages can be hidden and revealed plausibly.

## Variable Consistency

The Python implementation carefully maps to the theoretical concepts introduced in the PRIDE paper:

*   **`p`**: The large prime defining the finite field \( \mathbb{Z}_p \), used for all modular arithmetic.
*   **`n`**: The threshold for polynomial secret sharing (polynomial degree `n-1`), consistent with the paper's \(n\)-value.
*   **`K`**: The number of decoy messages generated per real encryption.
*   **`m`**: Represents a message (real or decoy) in the encryption/decryption process.
*   **`D`**: The FBE ciphertext, represented as a list of field elements (shares).
*   **`r`**: The key base, a list of random field elements used in the FBE construction for secret key generation.
*   **`sk`**: A secret key, represented as a tuple `(reconstruction_indices_D, chosen_indices_from_r)`. `reconstruction_indices_D` are the `n` x-coordinates whose corresponding y-values from `D` can reconstruct the polynomial, and `chosen_indices_from_r` are indices into the key base `r`.
*   **`params`**: A dictionary storing system parameters such as `p`, `n`, `K`, and `r_size`.
*   **`xs`, `ys`**: Lists of x-coordinates and y-coordinates, respectively, used in Lagrange interpolation.
*   **`coeffs`**: Coefficients \(a_1, \dots, a_{n-1}\) of the polynomial \(p(x)\) found by solving a linear system.
*   **`_encrypt_single`**: Maps to a simplified version of `FBE.Enc` incorporating polynomial sharing.
*   **`PRIDE_enc`**: Higher-level encryption function from the paper's `PRIDE.Enc` algorithm.
*   **`PRIDE_dec`**: Corresponds to `PRIDE.Dec` algorithm from the paper.
*   **`retrieval_mechanism`**: Corresponds to the paper's `M(D)` retrieval mechanism.

## Future Work (as per PRIDE paper)

*   **Adaptive Query Resilience:** Further research into mechanisms like dynamic threshold adjustment, sophisticated share redundancy, or verifiable random functions to strengthen PRIDE against adaptive attacks.
*   **Advanced Auditing Mechanisms:** Development of privacy-preserving auditing protocols leveraging Zero-Knowledge Proofs (ZKPs) and Verifiable Computation (VC) to address the "auditing paradox."
*   **Hybrid Privacy Models:** Exploration of integrating PRIDE with other privacy-preserving techniques such as Secure Multi-Party Computation (MPC), Homomorphic Encryption (HE), or Local Differential Privacy.
*   **Dynamic Data Handling:** Optimization of the underlying FBE scheme for efficient dynamic data management (robust share re-usage, incremental polynomial interpolation, verifiable data deletion).
*   **Side-Channel Mitigation:** Investigation and development of defenses against potential side-channel attacks targeting FBE implementation characteristics.
*   **Broader Applications:** Extending PRIDE's applicability to emerging domains like federated learning, IoT data aggregation, and decentralized identity systems.

---