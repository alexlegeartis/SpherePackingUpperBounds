# Independent audit: stronger dimension-five pd2 refinement

PASS: the fixed refinement gives density strictly below 0.510584101. All five verification stages were independently rerun against an isolated snapshot; no defect found. The published local-count certificates remain external mathematical inputs.

Source prefix: continuum_three_point/count_family/berncounts5pd2_R1.5_J1_p3_K12_N15_margin1e-07.
Isolated snapshot: continuum_three_point/count_family/audit_d5_pd2/candidate.json.
Parameters: dimension 5, R=3/2, J=1, radial feature degree 3, angular degree 12.

## Verification reruns

From workspace root, run continuum_three_point/count_family/certify_psd.py against the snapshot; run certify_count_pair_arb.py, certify_far_bernstein.py and certify_cyclic_bernstein.py from that directory against candidate_psd_exact.json. Run sphere_targets/certify_radial_dimension.py against candidate.json with argument 5.

Results, stored alongside snapshot:
- 13 exact rational positive-definite matrices; minimum reported LDL pivot 1.997697193311586e-8.
- Fourier Sturm: no positive roots, positive at zero.
- Physical-tail Sturm: no roots after conservative rational threshold below pi*R^2, negative at threshold.
- Count pair: 123 boxes, maximum depth 9, verified.
- Far Bernstein: 26211 boxes, 6971 sign-passed, 6135 excluded, maximum depth 45, verified.
- Cyclic Bernstein: 63 boxes, 28 sign-passed, 4 excluded, maximum depth 15, verified.

The reduced discovery PSD floor does not itself prove or weaken feasibility. The exact matrices are independently positive definite and are the same matrices used by all kernel sign certificates.

## Mathematical and domain review

The detailed dimension-five formula, support, endpoint and theorem analysis in continuum_three_point/count_family/audit_d5/AUDIT.md applies unchanged to this candidate. Verifier source and model parameters are unchanged; only rational coefficients/matrices/weights differ.

The normalized Gegenbauer recurrence is correct for S^4. The radial Gaussian-Laguerre basis has parameter 3/2 and Fourier eigenvalues (-1)^k. The exact rational origin values give f(0) and fhat(0), not floating solver normalizations. The function is Schwartz. Kernel positivity follows from positive-definite angular kernels and exact PSD coefficient matrices, with supported radial features and a separate degree-zero origin coordinate.

The cyclic numerator is (abc)^12 times the cyclic sum, with strictly positive multiplier on [1,3/2]^3. Every triple in this cube is a geometric triangle. Permutation exclusions use strict separation, retaining equality faces. Far triangles with exactly one side exceeding R leave a single kernel value; the far box covers all their angles and retains distance=R. Two long sides give zero. Pair subintervals include count cutoffs with the smaller right-hand allowance. Origin, r=1, r=R and one-sided support changes introduce no gap.

Exact rational initial Bernstein tensors followed by outward-rounded de Casteljau arithmetic enclose every checked box. The observed maximum depth 45 remains below the binary64 dyadic-endpoint limit. Far domain exclusion uses a rigorous maximum of squared distance over radial corners and the lowest cosine, with strict <R^2 exclusion.

Summing PSD local kernel inequalities and regrouping unordered triangles gives the same count-adjusted theorem:
density <= (pi^2/60)*(f(0)+K(0,0)+sum_i B_i lambda_i)/fhat(0).
The compactly supported kernel need not be continuous; finite local sums, periodic approximation and Schwartz radial decay suffice. There is no missing factor two on the count term.

## Neighbor-count assumptions

All four weights are nonnegative rational JSON decimals. Counts 45,46,47,48 are used at radii 1.003012,1.008096,1.013258,1.017445. Section 6/Table 2 of https://arxiv.org/html/1206.2608 supplies rigorous half-open counts at radii larger by 0.000001. Monotonicity gives the candidate's bounds. This audit checks correct use, not independent reproduction of the original published SDPs.

## Independent exact objective

fhat(0) = 500000000000006218545861/500000000000000000000000.
Unnormalized numerator = 15519895607095944475952251483/5000000000000000000000000000.
Normalized numerator = 15519895607095944475952251483/5000000000000062185458610000.

256-bit Arb density enclosure:
0.5105840999607987727021021850536941337926838178193477371147900754213375490648 +/- 5.44e-77.

Hence density <0.510584101.

Updated baseline: https://arxiv.org/html/2512.04829 (v2), Table 1, reports 0.5122645 in dimension five. Compared with this baseline and lower bound pi^2/(15 sqrt(2)):
- Precise certified enclosure closes 3.574795434244754% of the additive gap.
- Outward-rounded published bound 0.510584101 closes 3.574793223502243%.

These percentages compare against the displayed newer baseline, not an undisplayed exact value. The result is not the requested 20%-gap improvement.

Scope: source review, independent snapshot reruns and objective reconstruction, not clean-room reimplementation or formal verification. External local counts are relied upon as published results. No unrestricted literature-priority assertion is made.
