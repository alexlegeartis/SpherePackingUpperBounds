# Complete independent 16D count/profile density audit

## Verdict

**PASS: density in dimension 16 is strictly less than 0.024920988 for the frozen witness.**

The underlying expression is 0.0249209879146585407778922410323795054..., using the proved integer cap 7802. The conservative cap 7803 would instead give 0.0249210020077667745...<0.024921003. No change to the radial polynomial or kernel is needed when using the smaller proved count.

This is an independently checked research certificate, not external peer review or a formal proof-assistant verification. No numerical optimizer status or sampled constraint is a proof input.

## 1. Positive-width spherical count

The frozen candidate and factors have hashes b3225ad4cce819395ce3bd47613f15697af4646fb722ded1916020b5e52ba573 and 6860eca2c82fdb1b1ebb5d3a4108a216d46c900d6fac7bdeb0d880272e027954.

The local Gram blocks are exactly L^T L+lambda I, not L L^T. All local reweightings and global/nonorigin Schoenberg coefficients are nonnegative. The transverse homogeneous Gegenbauer factor is in dimension 15; latitude polynomials use dimension 16+2k. The separate origin and the degree-zero latitude feature have the stated cross and diagonal values. The nonorigin R kernel is PSD with zero origin row.

A fresh temporary-copy replay reconstructed all pair coefficients directly from these factors and checked the complete degree-32 pair interval by exact rational Sturm. A second reconstruction formed the entire cyclic polynomial from bilinear factors, checked permutation symmetry, and matched all 2673 coefficients exactly. A fresh exact integer replay then checked all 2155 proof events, covering the whole normalized cube with volume 1. It includes 641 direct sign leaves, 259 Gram-multiplier leaves, 107 infeasible-Gram leaves and 71 strict symmetry exclusions.

The Gram multiplier sign is correct: if G>=0, lambda>=0, and H+lambda G<0 on a box, then H<0 there. The proof does not discard Gram-boundary equality or sorted-boundary equality. Exact symmetry permits reduction to one ordered chamber. There is no reliance on the original contact SOS certificate or an unquantified perturbation of a contact-only bound.

The finite-sum spherical argument gives the real code bound 7802.215167128723216019597445... at cap 5003/10000. Integer cardinality therefore gives N<=7802.

## 2. Count-to-shell conversion

For a packing of minimum squared separation s=2, take all neighbors with squared radii in [2,tau), tau=10000/4997. Their normalized angular cosine is at most

    (u+v-2)/(2 sqrt(u v)).

The derivatives have the signs u-v+2 and v-u+2, which are positive since tau<4. Thus the maximum is attained at u=v=tau and is exactly 1-1/tau=5003/10000. The positive-width spherical cap gives at most 7802 neighbors in this shell. This is a pointwise cardinality bound, so integer flooring is legitimate. The radial pair proof uses the smaller, zero count allowance at the excluded outer boundary.

## 3. Kernel geometry and normalization

The scalar kernel uses dimension 16 Q with constant Schoenberg coefficient1, B=Q(1), and h(r)=exp(-phi(r^2/2)). Its exact entries are

    K00=C,
    K(0,x)=-(C/B)h(|x|),
    K(x,y)=(C/B^2)h(|x|)h(|y|)Q(cos(x,y)).

The degree-zero block is the outer product C(1,-1/B)(1,-1/B)^T. Higher blocks are nonnegative Schoenberg multiples. Therefore every local kernel sum is nonnegative. The diagonal is (C/B)h^2, so the pair correction is -(C/B)(2h-h^2), with no missing B or factor2.

The sibling geometry audit uses a separate closed-form Gegenbauer reconstruction, exact cap and root checks, a reconstructed binary cover, and analytic-gradient enclosures for every accepted rectangle. All 31,138 rectangles pass without further refinement. The all-scale argument remains valid at log/linear hinges by locally Lipschitz derivative hulls and the monotonicity of v phi'(v). The signed mixed-angle case and the complete infinite tail were checked explicitly.

## 4. Fourier, pair and rounding

The sibling radial audit independently reconstructed both the source Laguerre expansion and the direct monomial Fourier map. The epsilon repair and common normalization are exact. Rational Sturm proves the whole Fourier-positive ray and physical-negative ray. A new first-derivative interval proof covers all 59 finite pair intervals, all profile breakpoints and both count-threshold limits. The remaining infinite tail follows from the stronger physical-negative ray.

With the Fourier origin normalized to1, periodic summation gives

    rho <= P(0)+C+7802 lambda.

Multiplying by the volume pi^8/(256*8!) gives the reported density. This is not the unit-ball volume or the Barnes-Wall density. Rational Machin bounds independently prove outward rounding.

For clarity, adding the local PSD sums to the Poisson lower bound produces the upper expression N(P0+C) plus ordered pair allowances and twice each unordered-triangle cyclic sum. The cyclic terms are nonpositive, and the ordered shell degree is at most 7802 per center. Thus the count multiplier appears once, not twice. Gaussian outer decay makes all local periodic sums absolutely convergent. Periodic approximation extends the bound to arbitrary packings. No smoothness or Fourier transform of the isolated-origin kernel is required; only the radial Gaussian-polynomial f is used in Poisson summation.

## 5. Comparator and limits

On 2026-11-11 the primary arXiv HTML for version2, Table1, dimension 16 was checked directly:

https://arxiv.org/html/2512.04829v2#S5

It displays 0.02492121, whereas the older Cohn-et-al. value in the same row is0.02499441. We compare against the newer value. Our rounded bound is lower by 0.000000222, approximately 0.0008908% of that upper bound. The paper displays only eight decimal places; no unreported extra precision is invented. The gap exceeds ordinary last-digit rounding uncertainty by a large factor. The general live Cohn page failed to fetch during this check, so no claim of a new live-table update is made.

The gain is real but small. This audit makes no saturation or universal optimality claim, and no claim that the lower-density gap is substantially closed.

## 6. Safety checks and immutable inputs

The independent radial and geometry suites include ten negative controls. The unrepaired radial witness has a rigorously positive pair value and is rejected. Missing cover cells, invalid sorted exclusions, wrong count endpoint semantics, wrong repair sign, negative profile coefficients, overaggressive rounding and Python optimization are rejected.

The only review issue found was a missing no-optimization guard in one preexisting count replay entrypoint. Its owner added the guard; mathematical candidate, polynomial and proof-event bytes did not change. All newly supplied independent scripts already reject optimized Python.

No frozen producer mathematical inputs were changed during this audit.
