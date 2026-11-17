# Independent review of reconstructed 51D power-six candidate

**Rejected candidate. No density upper bound follows.**

The producer's exact edge checker found, and this audit independently reproduced without the producer's Jet or gate code, a positive cyclic kernel sum at squared side lengths `(1,1,67/52)` (cube `(0,0,5/8)`). The 512-bit lower endpoint is positive; raw value is about 7.19105388238. This is a genuine feasible triangle, not a parametrization singularity. Input SHA-256: `166963f7f6c8511a11ab9ab40ccce01f5cbc5edc12a32feb4557dba05a6192d8`.

## Scope of checks that passed

- All rational block principal minors are positive, and the chosen origin constant exceeds the exact Schur complement.
- Geometry and scalar kernel input bytes coincide.
- `H=(beta^6,beta^10,beta^14)` independently expands to the degree-28 scalar correction stored in the annular proof, with coefficients `2z_i/q` and `sum_k c_k B_k[i,j]/q^2`.
- The scalar Fourier-positive polynomial is byte-identical to the preceding published proof source.
- All 633 annular endpoints, widths, nonnegative heights, and cap-cost bindings are consistent. The conditional volume normalization is `1.3598001953322632e-7`, not a bound because the triangle premise fails.

## Parametrization review

For sorted squared sides `a=1+e x`, `b=a+e(1-x)y`, `c=b+e(1-x)(1-y)z`, `e=6/13`, all triples are feasible triangles because each length lies between 1 and `sqrt(19/13)<2`. After factoring the common nonnegative radial product, the cyclic expression has factors `(1-z)^6` and `((1-y)(1-z))^6`, with intermediate monomials in `beta^4`. The producer implements those factors correctly. At vanishing profile boundaries the original cyclic expression is zero, while the normalized sufficient condition extends continuously. The observed counterexample lies at nonzero profiles, so it cannot be dismissed as a boundary-extension issue.

The six Bernstein angular controls suffice for support-crossing triangles: the surviving angle lies in `[-1,1/2]`. The exact root-isolation checker passes these controls for this input. That does not imply the in-support three-side cyclic condition.

The standard positive-definite Gegenbauer addition theorem, kernel counting argument, Fourier/Poisson step, and periodic approximation remain mathematical premises. This review does not formalize them. No full rejected-geometry replay is needed: a certified counterexample already disproves its required sign.
