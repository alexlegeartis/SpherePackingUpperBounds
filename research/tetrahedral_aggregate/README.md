# Verified tetrahedral aggregate inequalities in dimension 22

This is a research-only, independently reproduced four-point proof package. It does not establish a new sphere-packing density bound, a packing construction, an optimal auxiliary constant, or literature priority.

For any finite unit-separated packing in R^22, let E count unordered edges with squared length at most q and T count unordered triangles with all three edges at most q. The package proves 3T <= BE for each of these two exact constants:

| q | Exact B | Strict upper | Prior stored affine certificate upper |
|---|---|---|---|
| 21/20 | 208821822550866169/62500000000000 | 3341.14916082 | 3846.43267004 |
| 11/10 | 1271745209351204667/250000000000000 | 5086.98083741 | 8042.25959581 |

The second constant is approximately 36.74687% smaller than the stored affine certificate. These comparisons are between valid witnesses, not proved optima. Exact incompatible endpoint slopes show that both new polynomials lie outside the entire affine-envelope family. NORMAL_HARMONIC_ENVELOPE_REDUCTION.md shows that adding unweighted higher face-normal harmonics cannot remove that obstruction in this shell-width range. Full joint tetrahedron geometry therefore supplies a genuinely larger feasible polynomial family.

## Replay

Python 3 and python-flint 0.9.0 are sufficient. No optimizer, network access, repository-relative external data, or previous packing certificate is needed.

    python -m pip install -r requirements.txt
    python verify_all.py

Do not use Python -O. The runner checks frozen-file hashes and exact agreement of independent witnesses, then runs the producer proofs, the two independent midpoint proofs, an alternative vertex-coordinate proof, the exact affine baseline signs, and endpoint separation. It uses a temporary copy, so frozen inputs and original audit logs are not overwritten. Expect a few minutes on one CPU. Every acceptance uses exact rationals or outward-rounded Arb intervals. Floating scores only choose the next subdivision.

## Contents

- PROOF.md: normalization, four-point summation, support and boundary conventions.
- SHARP_PROJECTION_CAP.md: independently reviewed far-angle geometry.
- NORMAL_HARMONIC_ENVELOPE_REDUCTION.md: exact all-degrees redundancy theorem.
- primary/: frozen rational certificates, producer interval implementations, exact baseline and separation checks.
- independent_q21_20/: explicit Gegenbauer reconstruction, independent midpoint and vertex-coordinate proofs, audit.
- independent_q11_10/: independent sequential monotone-face proof and audit.
- provenance/: discovery scripts/logs and the negative density activation report. Nothing here is a proof input.

The primary replay reconstructs exact coefficients from the two frozen decimal discovery rows and compares them with the independently frozen rational witnesses. The independent implementations read exact rational certificates directly. The primary affine baseline helper differs from its original source only in its default input filename, selecting the two local baseline rows.

## Density limitation

A bounded 22D density optimization using both improved constants chose effectively zero multipliers. Its numerical candidate remained noncompetitive and was not fully certified. The package must not be registered as a new density bound. A more general density kernel or weighting scheme would need a separate search and complete proof.
