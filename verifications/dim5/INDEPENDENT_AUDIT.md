# Independent dimension-five audit

PASS: the fixed candidate gives packing density strictly below 0.511534153, subject to the published local-neighbor bounds identified below. It also improves the newer 0.5122645 bound reported in arXiv:2512.04829v2, not merely the older 0.5126451306253027 benchmark. This is not a 20%-gap improvement.

## Fixed artifact and reruns

Source: continuum_three_point/count_family/berncounts5_R1.5_J1_p3_K12_N15_margin1e-07.json. An isolated copy is audit_d5/candidate.json; all rerun certificates alongside it point to that copy, not a moving producer file.

The portable package replays these checks with `python verify.py` against witness.json, reproducing the frozen rational matrices and asserting every result. The original workspace paths below identify provenance, not runtime dependencies.

Results:
- 13 exact rational positive-definite matrices; minimum reported LDL pivot 1.0000491026051282e-7.
- Fourier polynomial has no positive roots and is positive at zero.
- Physical polynomial has no roots past the conservative tail threshold and is negative there.
- Pair: 129 boxes, maximum depth 9, verified.
- Far: 1925 boxes (745 sign-passed, 218 excluded), depth 33, verified.
- Cyclic: 51 boxes (22 sign-passed, 4 excluded), depth 13, verified.

## Dimension-specific formulas

The normalized spherical polynomials satisfy P_0=1, P_1=t and (k+n-2)P_(k+1)=(2k+n-2)tP_k-kP_(k-1). For n=5 these are normalized Gegenbauer polynomials of parameter 3/2, positive definite on S^4. Both far and cyclic verifiers implement this recurrence correctly. The cyclic numerator expansion uses the coefficient of t^m divided by 2^m before expanding (a^2+b^2-c^2)^m, so it correctly accounts for cos(theta)=(a^2+b^2-c^2)/(2ab). Its multiplier (abc)^12 is positive.

The radial basis is e^(-pi*r^2)L_k^(3/2)(2*pi*r^2)/L_k^(3/2)(0). The implemented Laguerre recurrence and origin divisors are correct. Fourier eigenvalues are (-1)^k. Exact coefficients define a Schwartz function. The global Sturm certificates and Arb pair verifier use these same rational coefficients.

The rationalized matrices themselves define the kernel. All sign verifiers consume them, so they need not equal the floating solver matrices. The count weights are exact positive JSON decimal rationals.

## Continuum domain and theorem

The argument is the dimension-independent PSD cyclic packing-sum proof described in sphere_targets/audit_d4/AUDIT.md, with volume replaced by pi^2/60. Summing kernel PSD inequalities over centers leaves nonnegative total origin/diagonal contribution because each unordered triangle contributes twice its nonpositive cyclic sum. Adding the count-adjusted pair inequality and nonnegative-Fourier Poisson lower bound gives the numerator f(0)+K(0,0)+sum lambda_i B_i, divided by fhat(0). No additional factor of two belongs on the neighbor-count term.

The kernel is supported on {0} union [1,3/2], including annulus endpoints, with higher angular degrees vanishing at zero. This discontinuity is permitted by the direct finite-neighborhood sum proof and periodic approximation; no smooth-kernel assumption is needed.

All-short triangles are covered by the full [1,3/2]^3 cube modulo exact symmetry. Every such triple satisfies the triangle inequality. The symmetry exclusion is strict and retains equality faces. Exactly one long side leaves one kernel term, covered by the far certificate. For radii r,s<=R and distance>=R, cosine<=1/2: r^2+s^2-rs<=R^2 follows from separate convexity and corners. The far verifier's angular upper endpoint 1-1/(2R^2)=7/9 exceeds 1/2. Its distance exclusion uses strict <R^2, retaining the support boundary. Two long sides make all terms zero.

Pair subdivisions are closed. At count cutoffs the interval on the right uses the smaller allowance and includes the endpoint, proving the half-open count inequality. Beyond R the physical-tail certificate applies. At r=1 and R the pair certificate applies. Bernstein initial coefficients are exact rational; protected floating-point additions and halvings enclose subsequent de Casteljau coefficients. Observed maximum depths 33 and 13 preserve exact binary representation of domain endpoints.

## Published inputs and latest comparison

https://arxiv.org/html/1206.2608, Section 6/Table 2, provides dimension-five neighbor counts 45,46,47,48 for outer radii 1.003013,1.008097,1.013259,1.017446 and states rigorous verification. The candidate's four radii are each smaller by 0.000001, so the published half-open bounds apply. This audit verifies their use, not the original external SDP certificates.

https://arxiv.org/html/2512.04829v2, Table 1 (also available at https://arxiv.org/html/2512.04829), reports 5D upper bound 0.5122645. Therefore an improvement comparison must include that newer figure. A supplemental web search and the Cohn table https://cohn.mit.edu/sphere-packing/ did not identify a stronger finite 5D upper bound; that is not an exhaustive literature guarantee. We have not audited the newer paper's own numerical certificates.

## Exact objective

Independent exact rational reconstruction gives
normalization = 500000000000000874280347/500000000000000000000000,
normalized numerator = 3109754746473319928260607520391/1000000000000001748560694000000.

Multiplication by pi^2/60 using 256-bit Arb gives
0.5115341522016924552531670689836802262599793933614382953000898948389217703656 +/- 7.14e-77.

Hence density <0.511534153. Compared with 0.5122645 and the lower bound pi^2/(15 sqrt(2)), this closes approximately 1.55370382878685% of the gap. It cannot be advertised as achieving 20%.

Scope: verifier source review plus independent reruns and objective reconstruction, not a clean-room implementation of every verifier or a formally verified proof. No defect found in the candidate, continuum coverage, dimension adaptation, count aggregation or normalization.
