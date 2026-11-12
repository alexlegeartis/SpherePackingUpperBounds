# Reconstructed 16D angular refinement

**Strict upper bound: 0.024920883**, conditional only on replaying the included
linkage and the existing published base certificate as described below.

This is a new reconstruction from commit 043b119. It does not reproduce the
missing, stronger nested-shell package or its claimed 0.024899515 bound.

Place this package at `verifications/dim16_reconstruction_minB` alongside the
unchanged base `verifications/dim16`.

Run `python verify_all.py` from the repository checkout. This replays the
original base primary and independent chains, new primary and independently
implemented geometry, and both exact transfer calculations. The component commands are:

```
python verify_geometry.py
python verify_scaled_bound.py
```

Dependencies: python-flint. Search additionally uses NumPy/SciPy.

## Idea

Minimize the angular diagonal B=Q(1) while keeping the published radial
profile fixed and imposing the entire signed cyclic triangle inequality.
Exact rational coefficients and complete Arb interval covering establish
that the new angular kernel is admissible. Its nonnegative Schoenberg
coefficients prove positive semidefiniteness.

Replace C by C B_new/B_old. Then C/B is exactly unchanged, so the complete
published radial pair inequality, Fourier positivity, and count certificate
remain valid. Only the density objective's C term decreases. The independent
geometry and radial audit accompanies release integration.

The resulting density expression is
0.024920882308265284745..., strictly below 0.024920883.

Finite LP samples are discovery inputs only. `verify_geometry.py` covers
the full finite triangle domain and proves the infinite tail, while
`verify_scaled_bound.py` checks exact algebraic linkage and outward rounding.
The original radial/count replay is a required dependency, not replaced by
the geometry check. The supplied search output is not itself a proof.

The independent geometry audit checks all 17,637 accepted rectangles and the
35,925-node binary covering tree. Independent transfer rounding uses rational
Machin bounds for pi, rather than the primary Arb objective evaluation.
The stronger conditional radial research is deliberately excluded from this
release because its new spherical-count proof is not yet complete.
