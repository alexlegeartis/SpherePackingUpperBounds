# Independent computational audit of the four-dimensional witness

Result: PASS for a computer-assisted upper bound strictly below 0.635766242639, using the published local-neighbor bounds cited below. This does not attain the stated 20%-gap target. This audit reviewed the verifier implementations and reran them directly against an isolated snapshot; it did not independently reimplement every verifier or reproduce the published neighbor-count SDPs.

## Reproduction and exact witness

The frozen radial data and exact rational matrices are witness.json and witness_psd_exact.json. The audit reviewed verifier implementations and replayed them in a separate directory. Portable replay: python verify.py. Individual commands from this package directory:

python certify_psd.py witness.json
python certify_radial.py witness.json
python certify_count_pair_arb.py witness_psd_exact.json
python certify_far_bernstein.py witness_psd_exact.json
python certify_cyclic_bernstein.py witness_psd_exact.json

All pass. Exact rational LDL: 17 positive-definite matrices, smallest reported pivot approximately 9.92481593359753e-10. The rational matrices, rather than floating solver matrices, define the witness. Rationalization need not preserve the original solver values because all sign checks use the reconstructed matrices. Degree-31 radial coefficients and count weights mean their exact JSON decimal rationals.

Fourier polynomial: zero positive roots, positive at zero. Physical polynomial: zero roots after a rational threshold slightly below pi*R^2 and negative at that threshold. The Sturm algorithm rescales each remainder only by a positive factor; its endpoint sign-variation logic is appropriate here. The Gaussian-Laguerre radial function is Schwartz, so Poisson summation is justified. Its Fourier value at zero is exactly 249999999999726593313/250000000000000000000, not exactly one; the final quotient uses this exact normalization.

Pair certificate: 115 boxes, maximum subdivision depth 8. It evaluates the exact rational Chebyshev diagonal/cross polynomial together with the Gaussian-Laguerre radial function, using Arb second-order Taylor bounds. The derivative majorant is valid: the power-series coefficient absolute-value recurrence bounds each polynomial on its centered interval, and the Gaussian factor is bounded above by its value at the left radial endpoint. All count weights are nonnegative.

Far certificate: 34009 boxes, 8439 sign-passed and 8566 excluded, maximum depth 29. Cyclic certificate: 83 boxes, 38 sign-passed and 4 excluded, maximum depth 11. Both terminate with no unprocessed boxes. Their exact rational Bernstein initialization and nextafter-protected additions/halvings enclose de Casteljau coefficients; at these depths all dyadic domain endpoints remain exactly representable in binary64.

## Domain and boundary coverage

Kernel features are defined to vanish off {0} union [1,R], with R=333/200=1.665, and use the closed annulus including its endpoints. Higher angular-degree features vanish at the origin. The angular basis U_k(t)/(k+1) is positive definite on S^3, normalized to one on the diagonal. Positive-semidefinite coefficient matrices therefore give a positive-semidefinite supported kernel, even though it is discontinuous at the support boundary.

If all triangle sides are in [1,R], the cyclic verifier covers their whole cube modulo permutation symmetry. Since R<2, every triple in that cube is a genuine nondegenerate triangle. The expanded polynomial is (abc)^16 times the sum over its three vertices; the multiplier is strictly positive. Expansion coefficients follow U_k(t)=sum_j (-1)^j binom(k-j,j)(2t)^(k-2j). Symmetric Q matrices make the resulting cyclic expression symmetric under all side permutations. Domain exclusion uses strict separation from the ordered region, so equality boundaries are retained.

If precisely one side exceeds R, only the kernel term using the other two sides remains. The far verifier proves that term nonpositive for distance >=R, including the boundary. Its angular box extends from -1 to 1-1/(2R^2), which covers all relevant angles: for 1<=r,s<=R and distance>=R, cos(theta)<=1/2. Indeed r^2+s^2-rs<=R^2 by separate convexity and checking the four corners. The actual upper endpoint used is larger than 1/2. Far domain exclusion bounds squared distance above by evaluating radial corners at the smallest angle cosine; separate convexity makes this sound. Strict <R^2 is necessary and is used, retaining the boundary. Two or three sides exceeding R give zero cyclic contribution.

At a count radius, the right-hand subinterval checks the smaller allowance on its closed left endpoint. This proves the required half-open condition. At r=1 and r=R the closed pair certificate applies. Beyond R the kernel is zero and the global physical-tail sign applies. Thus no one-sided support or count-cutoff gap remains.

## Packing theorem and counts

For a periodic packing, sum the PSD kernel over the finite displacement neighborhood around each center, including zero. Separate origin and diagonal terms. The remaining off-diagonal terms regroup into twice each triangle's cyclic sum, which is nonpositive. Consequently the total origin-and-diagonal contribution is nonnegative. Combining this with the count-adjusted pair inequality and the nonnegative-Fourier Poisson inequality gives

density <= (pi^2/32) * (f(0)+K(0,0)+sum_i B_i lambda_i)/fhat(0).

One may pass to arbitrary packings by large finite truncations and periodic replication separated by a fixed-width gap. Compact support bounds all kernel boundary costs. The Schwartz radial function permits the usual Poisson/periodic approximation argument. This direct proof does not require continuity of the supported kernel.

Published input: David de Laat, Fernando Mario de Oliveira Filho and Frank Vallentin, Upper bounds for packings of spheres of several radii, https://arxiv.org/html/1206.2608, Section 6 and Table 2. Table 2 gives dimension-four counts 24,25,26,27 at outer radii 1.008097,1.017446,1.025978,1.036951. Section 6 defines the relevant neighbors at distances in [1,1+epsilon) and states rigorous verification. The candidate uses each radius smaller by 0.000001, so these bounds apply monotonically. This audit relies on that published result, not on a fresh reproduction of its SDP certificates.

## Objective and comparison

At 256-bit Arb precision the exact reconstructed objective is enclosed by
0.6357662426387588618572394093758579880966147479317909183832447245178820105673 +/- 5.00e-77.

Thus the outward-rounded decimal 0.635766242639 is a valid concise upper bound. The expression uses rational f(0), rational Q_0[0,0], rational count contributions and exact rational fhat(0), multiplied by an Arb enclosure of pi^2/32. It is strictly below the supplied benchmark 0.6361073321551329, and strictly above the target 0.6322559207377233. Relative to that benchmark and lower bound pi^2/16, gap closure is approximately 1.771244250002509%. This audit establishes the comparison to the supplied benchmark; it does not independently survey whether that benchmark is the strongest published result today.

The package freezes the audited snapshot and its exact rational matrices. Later optimization output must not be substituted without rerunning every check. The repository reports the slightly more conservative decimal 0.635766243 in its comparisons.
