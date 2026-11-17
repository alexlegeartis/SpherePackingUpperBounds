# Verified lattice-only refinements,36–48 dimensions

**Method:** Dimension-adapted positive-definite angular kernels with universal rank-two lattice sums, full collinear-line domination, and conditioned mod2 class profiles. Dimension36 additionally uses a nonnegative quadratic Fourier multiplier and reoptimized class-cut weights.

These are lattice-only upper bounds. They do not replace unrestricted-packing upper bounds or establish optimality. [Proof](PROOF.md) · [Registered conservative bounds](REGISTERED_RESULTS.json).

Install requirements.txt and run `python run_all.py`. To replay one dimension use `python run_all.py --dimension 40`. Assertions must remain enabled. The runner verifies frozen input hashes, materializes hash-checked compressed JSON inputs in temporary storage, and reruns the relevant exact/interval sign gates, source radial/profile dependencies and final rounding. A full replay may take tens of minutes. No optimization solver or network is needed.

The final46D upper is7.858622654e−7 using a nonnegative quadratic Fourier multiplier;48D is3.527782891e−7 using a quadratic multiplier plus an additive Fourier sum of squares. Both retain the proved scalar rank-two kernel. No harmonic-dependent matrix candidate is included.

## Audit scope
The angular polynomial identities, PSD coefficients and dimension restrictions are exact rational checks. Compact domains use outward interval arithmetic; the radial replay uses Arb. Source class-cap and weighted-profile theorems are checked algebraically. Independent audits reran the scalar gates at60–100 decimal digits and radial checks at1024 bits; the original predecessor proof hypotheses and public mathematical inputs are stated in PROOF.md. A shared arithmetic backend is not claimed as an independent implementation.

The mod2 framework extends [Zubrilina(2026)](https://arxiv.org/abs/2610.35830); the Fourier/positive-definite method uses [Cohn–de Laat–Salmon](https://arxiv.org/abs/2206.15373). The new rank-two summation and collinear transfer are proved explicitly here. Each registered row compares to the preceding SpherePackings commit20b638f5d5b0bf00a2a4ce8d70c4e062517a502e; all original packages remain unchanged.
