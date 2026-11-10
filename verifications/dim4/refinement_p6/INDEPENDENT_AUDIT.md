# Independent audit: dimension-four p6 refinement

PASS: fixed audited_snapshot.json proves density strictly below 0.633430327, relying on the published neighbor-count inputs identified below. All five verification stages were independently replayed against an isolated snapshot. No mathematical or continuum-coverage defect found.

Snapshot: audited_snapshot.json.
SHA256: 149aeac1edf877ca937397495a58169d2f19c4770781050c8581420ae3b1308f.
Snapshot bytes match the producer source; reconstructed exact matrices match the producer matrices.
Parameters R=333/200, J=1, p=6, K=16.

## Verification replay

From this package directory, run python verify.py. Component commands:
python certify_psd.py witness.json
python certify_radial.py witness.json
python certify_count_pair_arb.py witness_psd_exact.json
python certify_far_bernstein.py witness_psd_exact.json
python certify_cyclic_bernstein.py witness_psd_exact.json

Results stored beside snapshot:
- 17 exact positive-definite rational matrices, minimum reported LDL pivot 8.873087728341943e-10.
- Fourier polynomial has zero positive roots and is positive at zero.
- Physical polynomial has zero roots after a conservative rational threshold below pi*R^2 and is negative at that threshold.
- Count-adjusted pair: 215 boxes, maximum depth 8, verified.
- Far: 17723 boxes, 4973 sign-passed, 3889 excluded, maximum depth 37, verified; angular_top=1/2.
- Cyclic: 135 boxes, 55 sign-passed, 13 excluded, maximum depth 13, verified.

The producer's early message said 139 cyclic boxes. Its finalized report and this independent rerun both say 135. This is a reporting discrepancy, not a feasibility failure.

## New angular cap: proved safe

For the far region let 1<=r,s<=R and d>=R. The law of cosines gives
cos(theta) <= (r^2+s^2-R^2)/(2rs).
To prove the right side <=1/2 it suffices to prove r^2+s^2-rs<=R^2. The left side is convex in each variable separately; its maximum on [1,R]^2 is attained at a corner. Corner values are 1, R^2, and 1+R^2-R (twice). Since R>=1, all are <=R^2. Thus every relevant far triangle has cosine<=1/2, with equality included.

The verifier change sets the exact rational angular cap to 1/2 rather than the previously broader cap. It retains the endpoint. Consequently this change discards no required configurations and makes certification easier without weakening its conclusion. The lower cosine endpoint remains -1.

## Remaining continuum and theorem checks

The detailed formula/theorem review in ../INDEPENDENT_AUDIT.md applies unchanged with radial degree six. In particular the radial function is a Schwartz Gaussian-Laguerre expansion with Fourier eigenvalues (-1)^k; all sign checks use exact JSON decimal coefficients. The rationalized coefficient matrices define the kernel, and all sign checks use these same matrices. Normalized U_k(t)/(k+1) angular kernels are positive definite on S^3. Supported radial features and the separate origin coordinate preserve kernel PSD, despite cutoff discontinuities.

Cyclic Bernstein checks cover [1,R]^3 modulo permutation symmetry, with positive multiplier (abc)^16. All such triples satisfy triangle inequalities because R<2. Symmetry exclusions are strict and retain equality boundaries. Exactly one side >R leaves one kernel term, covered by the far certificate; two long sides give zero. The distance exclusion is strict (<R^2) and therefore retains the d=R boundary. At R itself the full cyclic condition still applies and is certified.

The pair checker covers closed subintervals and assigns the smaller allowance at every count cutoff from its right-hand interval. r=1 and r=R are included, while the global physical-tail sign applies beyond R. There is no omitted one-sided endpoint condition.

Exact rational Bernstein initialization and outward-rounded float64 de Casteljau arithmetic enclose all coefficients. Actual maximum depth 37 is below the dyadic-endpoint precision limit. All stacks exhaust with no unresolved boxes.

The packing theorem follows by summing finite local PSD inequalities, regrouping each unordered triangle into twice its nonpositive cyclic sum, then combining the nonnegative origin/diagonal contribution with count-adjusted pair inequalities and Fourier Poisson summation. The resulting objective is (pi^2/32)*(f(0)+K(0,0)+sum B_i lambda_i)/fhat(0). Periodic approximation with fixed-width gaps and Schwartz decay extends the bound to arbitrary packings. No smoothness of the supported kernel is needed.

## External count input

The four nonnegative weights multiply counts 24,25,26,27 at radii 1.008096,1.017445,1.025977,1.036950. Section 6/Table 2 of https://arxiv.org/html/1206.2608 supplies rigorous half-open neighbor bounds at radii larger by 0.000001. Monotonicity therefore establishes the exact count assumptions used. This audit verifies their applicability, not the original published count SDP calculations.

## Exact objective and current comparison

Independent rational reconstruction gives:
fhat(0) = 399999999999996070879/400000000000000000000.
Unnormalized numerator = 8215028529461860492858846879/4000000000000000000000000000.
Normalized numerator = 8215028529461860492858846879/3999999999999960708790000000.

256-bit Arb density enclosure:
0.6334303260113454414754145742596239771574365762937475831059441425951904455188 +/- 8.16e-77.

Thus density <0.633430327.

Updated 4D comparison: arXiv:2512.04829v2 Table 1, https://arxiv.org/html/2512.04829, reports 0.63610277. Relative to the additive gap between that displayed baseline and pi^2/16:
- Precise enclosure closes 13.88102683888734%.
- Outward-rounded 0.633430327 closes 13.88102170368506%.

The result is still below the requested 20%-gap target. These calculations use the published displayed baseline, not its undisplayed exact value.

Scope: verifier source review, independent snapshot replay and exact objective reconstruction. Not a clean-room implementation, formal proof-assistant verification or independent reproduction of external neighbor-count SDPs. No unrestricted literature-priority assertion follows from this audit alone.

The portable witness.json retains exactly the mathematical fields of audited_snapshot.json; the runner asserts equality.
