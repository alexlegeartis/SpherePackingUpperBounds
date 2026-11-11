# Independent audit of the 6D cyclic/count trial

The final frozen 6D candidate passed a fresh complete producer replay and independent algebra checks. Earlier failed snapshots described below remain invalid; only the final certificate named here is endorsed. This is a source-reviewed computational certificate, not a formal proof-assistant verification or an independent implementation of every continuum checker.

## Analytic checks

For each center, sum the PSD kernel over its finite supported neighbor set together with the origin. Summing this nonnegative expression over centers produces one origin diagonal per center, one `2K(x,0)+K(x,x)` per ordered neighbor pair, and twice the three-vertex cyclic sum per unordered triangle. Adding the radial pair sum, applying the pair allowance and cyclic signs, and using Poisson summation gives precisely the stated numerator divided by the Fourier value at zero. Consequently the published average-degree costs have no additional factor of two. The radial function plus both Gaussian guards is Schwartz; the compactly supported kernel need not be smooth because Poisson summation is applied only to the radial function.

The normalized Gegenbauer recurrence used in both verifiers is correct in dimension n. Nonnegative PSD blocks imply a PSD kernel by the addition theorem; zeroing the feature map outside its radial support and isolating the origin preserve this property. A triangle with two edges outside support contributes zero at every vertex. A triangle with exactly one long edge has exactly one potentially nonzero vertex contribution. At the support boundary the compact checker includes all three terms, while the far checker imposes the stronger single-term condition as well.

For 1<R<2 every side triple in [1,R]^3 satisfies the triangle inequality. The cyclic numerator construction uses the correct positive denominator (abc)^K. Its three permutations are correct because each matrix is symmetric. The far angle interval includes all needed geometries: when the third side is at least R and both incident radii are at most R, their cosine is at most 1/2, which is no greater than the implemented upper endpoint 1-1/(2R^2).

Both reciprocal Gaussian transform formulas and origin costs are correct. Dividing the original Fourier sign inequality by exp(-Y) bounds its potentially negative slow guard below by -eta*2^n*(2n-c), as implemented. All omitted physical guard terms are nonpositive on the packing domain. Exact origin Schur saturation changes no off-origin constraints. Its correction must use the guarded normalization, as implemented.

The chosen reciprocal scales require n/8<c<pi. Thus this particular repair does not extend to n>=26. The 6D implementation is valid on that issue.

## Published dependency

Independently checked https://arxiv.org/html/1206.2608, Section 1.4, Section 6 and Table 2. Epsilon-tangency is defined by distances [1,1+epsilon), and the required dimension-six epsilon/cost pairs are (.002006,79), (.004024,80), (.006054,81), (.008097,82). The paper explicitly describes its bounds as rigorous. Decreasing each threshold by .000001 preserves the bound by monotonicity. These are externally published mathematical dependencies, not independently rerun spherical-code certificates.

## Implementation checks and hardening

The pair verifier encloses a complete second-order Taylor remainder and splits at each threshold. The right-hand subinterval includes each threshold and uses the smaller allowance. Its polynomial derivative and chain-rule factors are correct. The physical tail begins at R^2*pi_lower, a conservative earlier point than the true support endpoint.

The exact rational tensor-Bernstein initialization and outward-rounded float64 de Casteljau steps are conservative. NaN or infinity cannot accidentally pass the strict-negative test. Original geometry/symmetry boxes used floating midpoint coordinates with depth limits beyond binary64's mantissa capacity. This was reported; the parent hardened both verifiers to Fraction endpoints. I inspected the updated code: subdivision midpoints and symmetry comparisons are exact, and far geometric endpoints enter Arb through exact rational conversion. This resolves that potential coverage concern. The previously observed compact proof depth was only 16, so the issue was not evidence of an invalid earlier compact result.

## Independent arithmetic executable

`audit_independent.py EXPORTED_DIRECTORY` constructs Laguerre polynomials by explicit binomial coefficients and verifies the complete Fourier polynomial independently using repeated radial Laplacians. It independently checks the origin Schur complement, guarded normalization, and full saturated objective. It also constructs exact tensor-Bernstein coefficients of the antipodal kernel on an 8-by-8 partition of the complete radial square. This last test is sufficient when it passes, but failure of a positive Bernstein coefficient alone does not prove a positive function value.

On `witness6_full_mid`, Fourier and origin arithmetic pass. The antipodal sufficient test fails; the full far verifier also fails, as separately known. The objective approximately 0.410176576 is therefore not a density upper bound.

On `witness6_antipodal_mid`, all independent algebra checks pass, and the full antipodal boundary is certified nonpositive by exact Bernstein coefficients. The largest coefficient is exactly -1601176835260239/1600000000000000000000000 (approximately -1.000735522e-9). The independently calculated objective is 0.4101726628679345567106601023549950325..., still provisional because a full far-domain verification has not passed. Results are stored beside the witness as `independent_algebra_audit.json`.


## Final certificate audit: PASS

Frozen input: `candidate_certificate_v2_6.json`.
SHA-256: `797128e7ada65976ce25926ee84f626a9607208013bfb420243215df2f871606`.

I ran a fresh full producer replay from working directory `/tmp`, using absolute script/input paths and fresh output directory `certificate6_v2_independent_replay`. This exercised local module import portability and all three verifiers' relative `witness.json` resolution independently of the research working directory. All exact PSD, Gaussian sign repair, count inputs, pair, far and compact cyclic checks passed. The output log is `certificate6_v2_independent_replay.log`.

- Pair continuum: 197 checked boxes, maximum depth 11.
- Full one-long-edge continuum: 3,097 checked boxes, 1,135 sign-certified leaves, 414 geometrically excluded leaves, maximum depth 42.
- Full compact cyclic continuum: 153 checked boxes, 70 sign-certified leaves, 7 permutation-excluded leaves, maximum depth 14.
- Independent explicit-Laguerre/radial-Laplacian Fourier polynomial identity: PASS.
- Independent origin Schur complement, Gaussian origin values, and corrected normalization: PASS.
- Independent complete antipodal boundary: PASS, exact largest coefficient -95091611845759/100000000000000000000000.

The full fresh replay verified the strict upper decimal **0.410189987**, above the outward Arb enclosure centered at

0.4101899868311988677407604362994285017515171298713618040883764151908393153314

with reported radius 4.46e-77. The density normalization is the volume of the six-dimensional radius-one-half ball, pi^3/384, multiplying the rational saturated objective. The origin correction is divided by the repaired Fourier value, not the original one.

The six-dimensional frozen witness has R=3/2, p=4 and angular degree K=12. Exact published count thresholds are 1.002005, 1.004023, 1.006053, 1.008096, with weights nonnegative and bounds 79, 80, 81, 82. Their use is cumulative, with half-open upper endpoints, as required. No domain remains inferred from finite samples: all compact pairs, compact triangles and surviving one-long-edge configurations are covered, and both unbounded radial signs are checked separately.

Source review confirms the hardened Fraction boxes, conservative Arb geometric exclusion, and outward-rounded float64 Bernstein subdivision remain in the rerun. Assertions are essential to these programs. Explicit `python -O` refusal was tested successfully for the exporter, replay driver, all three continuum verifier entrypoints, and Gaussian guard helper. The independent audit executable now likewise refuses optimized mode. Details of the five main entrypoint tests are in `certificate6_v2_optimized_mode_tests.json`.

Thus the bound is below the comparison value 0.4102905 reported in Table 1 of https://arxiv.org/html/2512.04829v2. This audit does not claim a worldwide-priority result, exhaustive current-literature dominance, or independent reproduction of the external published spherical-code count certificates. It establishes successful source-reviewed replay of the supplied 6D certificate conditional on those valid published count inputs and the analytic theorem reviewed above.
