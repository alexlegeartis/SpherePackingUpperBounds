# Independent47D restriction and origin-Gram audit

The exact algebra verifier independently reconstructs the claim
Delta_47 <0.000000616879891.
The density enclosure is centered at6.16879890934120654755443575858423096817928440380657630319871654321838145113197674742058405146926974592292609929612e-7. All three copied source48 continuum verifiers passed freshly: angular exact checks; radial112-coefficient Fourier identity, root isolation and33 coupling intervals; and11,932 compact cyclic rectangles plus analytic tails. Their logs accompany this review.

## Restriction and Fourier normalization

Let f48(x)=P(pi|x|²)exp(−pi|x|²) use the unnormalized degree111 source polynomial. Set f47(x)=f48(x,0). Since f48 is Schwartz, its47D Fourier transform equals the integral of fhat48(xi,t) over the missing real frequency coordinate. Source Fourier positivity thus implies47D Fourier positivity without a new root-isolation claim. All integrals converge absolutely.

The independent verifier computes the full47D Fourier polynomial directly from Gaussian monomials, and separately computes the marginal coefficients by
h47_j = sum_(k>=j) h48_k binom(k,j)(1/2)_(k−j).
It checks equality of all112 coefficients, independently checks the48D monomial transform, and computes the positive normalization N47=sum p_k(47/2)_k exactly. Using the old48D Fourier origin would be incorrect; the objective here uses N47.

## Dimension-specific angular PSD and origin change

Q is the same exact degree19 source angular polynomial and B=Q(1). Explicit rising-factorial formulas, rather than the producer's Gegenbauer recurrence, construct normalized47D Gegenbauer polynomials. Exact triangular conversion yields Q=sum b_l P_l^(47), with every b_l>=0, sum b_l=B, and b0>1. An independent spherical moment identity checks b0=sum Q_(2j)(1/2)_j/(47/2)_j.

Set raw M=neighbor_cost/B*(1000001/1000000), kappa17/8, r²6, and h(x)=exp(−kappa*pi*(|x|²−6)). The non-origin kernel at allowed separations remains (M/B)h(x)h(y)Q(angle). Its degree-zero neighbor coefficient is M*b0/B. The origin/cross Gram matrix is
[[M*B/b0, −M], [−M, M*b0/B]].
Its positive diagonal entries and zero determinant are checked exactly. Therefore it is positive semidefinite. Higher harmonic terms are PSD by the47D addition theorem. The kernel has K(0,0)=MB/b0, K(x,0)=−Mh(x), and K(x,x)=Mh(x)² for |x|>=r. The latter two values and every non-origin triangle value are unchanged from the raw source.

A smooth global extension avoids an isolated-origin artifact: choose a smooth radial cutoff chi that is zero near0 and equals1 for |x|>=r. Use degree-zero feature sqrt(b0)*chi*h−B/sqrt(b0)*(1−chi), with overall factor M/B. Use chi*h times the positive-degree harmonic features. These are smooth, rapidly decreasing, and give exactly the required values. No new packing assumption is introduced.

## Inherited global inequalities

The source radial verifier proves f48<=M(2h−h²) at all radii>=sqrt6. Restriction preserves this inequality pointwise. The kernel change preserves cross and diagonal terms, so the coupled pair condition is inherited unchanged.

The cyclic expression depends only on three side lengths and the same angular polynomial, not on ambient dimension. The source compact interval proof and analytic tails therefore apply unchanged in47D. Its scale argument uses Q's single sign change and the ordering of the two secondary angles, not a48D-only geometric condition. In normalized variables, alpha=(51/4)*pi>40 is independently checked with rational pi bounds. For t<=0 both exponential gaps are at least1. For u>=5 a positive first secondary term forces its gap>4.001 and the second gap is at least24, yielding the source rational exponential bound. The compact domain covers the remaining sorted triangles.

The kernel is rapidly decreasing, so periodic pair/triple sums converge absolutely. The same finite PSD/cyclic argument and periodic Poisson summation bound number density by (P(0)+MB/b0)/N47. Periodic approximation extends to arbitrary packings. This is not a lattice-only result and not a claim that the kernel satisfies the literal pointwise-negative premise of Theorem1.1.

## Exact final rounding

The47D ball radius is sqrt6/2. Its squared volume is4*pi^46*6^47*(24!)²/(48!)². The independent verifier uses rational alternating-series Machin bounds with80 terms and proves the squared final density below(0.000000616879891)², separately from Arb display rounding.

This is a reproducible computer-assisted audit, not external peer review or proof-assistant formalization. A fresh search retrieved the indexed47D row at https://cohn.mit.edu/sphere-packing/: upper0.0000006490757338, lower0.0000000226696032…, center density2^(-118)*47^(47/2). The index is four months old and the live page returns502, so this is a verified pinned comparison, not an exhaustive current literature guarantee. No stronger finite47D general upper bound surfaced in targeted searches. No worldwide-priority claim follows solely from this calculation.
