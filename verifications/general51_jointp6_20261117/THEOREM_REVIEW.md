# Independent review: new joint 51D power-six certificate

Full geometry passed: all 53,227 leaves of the 106,453-node partition were replayed independently at 512 bits. Together with the exact algebra, support-crossing angular, full radial and infinite-tail checks, this certifies unrestricted 51D packing density < 1.37012082e-7.

## Scope and counting argument

This is an unrestricted sphere-packing argument, not a lattice or antipodal argument. For nonzero displacement vectors the compact PSD kernel has three features `H=(beta^6,beta^10,beta^14)`, `beta=max((19/13-r^2)/(6/13),0)`, and thirty normalized 51-dimensional Gegenbauer harmonic blocks. Every rational block is positive definite. The exceptional origin value and origin-to-nonzero feature vector satisfy the exact degree-zero Schur inequality. The Gegenbauer addition theorem therefore gives an origin-extended PSD kernel. This remains valid despite the kernel's separate definition at the origin: continuity there is not needed for finite PSD sums.

For three packing centers, the cyclic sum of the three vertex kernels must be nonpositive. For all three squared sides in `[1,19/13]`, sorted sides are parameterized by the full unit cube. The independently replayed sufficient expression factors out a nonnegative common radial product and uses the factors `(1-z)^6` and `((1-y)(1-z))^6`, with residual monomials in `beta^4`. Every such side triple is a feasible triangle since `sqrt(19/13)<2`; covering the full cube is sufficient. Boundary extension is conservative and cannot hide an interior positive value.

If exactly one side exceeds the support cutoff, the surviving angle cosine is at most `1/2`. All six exact tensor-Bernstein controls of the symmetric bidegree-(2,2) angular profile are strictly negative throughout `[-1,1/2]`, by derivative root isolation. The power-to-Bernstein map is `[[1,0,0],[1,1/2,0],[1,1,1]]`. If two sides exceed support, every cyclic term vanishes.

Sum the kernel PSD inequality over centers of a periodic packing. Distinct triples contribute twice their nonpositive cyclic sums. Consequently the averaged diagonal-plus-origin correction is at least `-s`. Add this correction to the scalar Schwartz function `exp(-X)P(X)`, `X=pi*r^2`, minimum separation `X=22`. Its positive residual is bounded annulus by annulus; all positive-height annuli use individually checked exact spherical-code LP polynomials. The radial widths satisfy the geometric condition for their angle bounds. Poisson summation then bounds density by

`(22/4)^(51/2)/Gamma(53/2) * (P(0)+s+budget)`.

Periodic approximation extends the conclusion to arbitrary sphere packings. Standard Gegenbauer positive definiteness, the spherical-code LP theorem, Poisson summation and periodic approximation are mathematical inputs, not formalized theorems in this package.

## Independent checks

- Rational principal minors and origin Schur complement, scalar/geometry input byte identity.
- Explicit Laguerre coefficient reconstruction of the new input Fourier polynomial. The new scalar is not the previous published scalar.
- Exact positive factorization of the repaired Fourier polynomial and an independent closed-coefficient Fourier involution, without relying on the producer's recurrence implementation.
- Direct expansion of the compact scalar correction, of degree 28, agrees with the annular certificate.
- Every coefficient of `P(X+180)` is nonpositive and the leading coefficient is negative, supplying a direct exact infinite-tail proof independent of real-root classification.
- Full fresh replay of all 633 radial annuli. The polynomial/exponential interval algorithm is shared with the producer, reviewed and rerun, not independently reimplemented. The regenerated annular digest is identical.
- All 11 positive-height annuli have exact polynomial code bounds checked independently using a different Gegenbauer expansion path. Independent 512-bit ball-volume normalization is strictly below the reported decimal.
- Complete geometric partition was replayed at 512 bits using all nine matrix entries and the generic `arb_poly` Taylor-shift Jet, rather than the producer's six symmetry-grouped entries and fast shift. Prefix audit output is never a completion claim.

The previously rejected `producer_p6` is a different candidate; its positive edge value is not present evidence for or against this candidate. Frozen hashes identify the exact new inputs.
