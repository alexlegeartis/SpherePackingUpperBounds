# Independent audit: dimension-five p4 refinement

PASS: fixed continuum_three_point/count_family/berncounts5p4_R1.5_J1_p4_K12_N15_margin1e-07.json yields density strictly below 0.510522525. All verification stages were replayed against an isolated snapshot. No defect found. Published local-count certificates remain external inputs.

Snapshot: continuum_three_point/count_family/audit_d5_p4/candidate.json.
SHA256: 8c423636bebe75aed179604a1ad65c5555dc85d5c8efb9d13faaff002957e1ff.
Parameters: dimension 5, R=3/2, J=1, radial feature degree 4, angular degree 12.

## Replay

From workspace root run the count_family scripts certify_psd.py on candidate.json, then certify_count_pair_arb.py, certify_far_bernstein.py and certify_cyclic_bernstein.py on candidate_psd_exact.json. Run sphere_targets/certify_radial_dimension.py on candidate.json with argument 5. Results are saved beside the snapshot.

- 13 exact rational positive-definite matrices; minimum reported LDL pivot 1.9751853601993075e-8.
- Global Fourier polynomial: zero positive roots, positive at zero.
- Physical tail: zero roots beyond conservative threshold below pi*R^2, negative at threshold.
- Count-pair: 127 boxes, maximum depth 8, verified.
- Far: 5373 boxes, 1523 sign-passed, 1164 excluded, depth 33, verified.
- Cyclic: 85 boxes, 39 sign-passed, 4 excluded, depth 11, verified.

## Formula, continuum and theorem checks

The detailed audit in continuum_three_point/count_family/audit_d5/AUDIT.md applies with feature degree four. The verifier source and dimension conventions are unchanged. In particular normalized Gegenbauer polynomials with parameter 3/2 are positive definite on S^4; the radial Gaussian-Laguerre basis also has parameter 3/2, with Fourier eigenvalues (-1)^k and correctly divided origin values. The rational matrices, not approximate solver matrices, define the kernel and are used throughout certification. The radial function is Schwartz.

The cyclic numerator is (abc)^12 times the cyclic sum, with positive multiplier on [1,3/2]^3. Every cube triple is a geometric triangle. Permutation exclusions are strict and retain equality boundaries. Exactly one side >R leaves one far kernel term; two long sides leave zero. Far angular cap 1-1/(2R^2)=7/9 covers the necessary cosines, which are <=1/2 by separate convexity of r^2+s^2-rs and corner evaluation. Distance exclusion uses strict <R^2 and retains the support boundary.

Pair subdivisions include their endpoints; the interval to the right of each count cutoff uses the smaller allowance at that endpoint. r=1 and r=R are included, and the global physical-tail sign applies beyond R. There is no uncovered boundary or one-sided condition. Exact Bernstein initialization plus outward-rounded de Casteljau arithmetic encloses every polynomial box. Observed depths 33/11 remain within exact binary64 dyadic-endpoint representation.

The PSD cyclic sum argument gives density <= (pi^2/60)*(f(0)+K(0,0)+sum B_i lambda_i)/fhat(0), with each unordered triangle contributing twice its nonpositive cyclic sum. The count term has no extra factor two. Kernel discontinuity at its compact support does not invalidate the finite local-sum proof; periodic approximation and Schwartz radial decay yield the general-packing result.

## Count assumptions

The four count weights are nonnegative. Counts 45,46,47,48 are used at radii 1.003012,1.008096,1.013258,1.017445. Section 6/Table 2 of https://arxiv.org/html/1206.2608 gives rigorous half-open counts at radii larger by 0.000001. Monotonicity makes those published bounds applicable. This audit verifies usage, not the original external SDPs.

## Independent exact objective and updated comparison

fhat(0) = 500000000000013867843581/500000000000000000000000.
Unnormalized numerator = 19397529933890261448094669/6250000000000000000000000.
Normalized numerator = 38795059867780522896189338/12500000000000346696089525.

256-bit Arb density enclosure:
0.5105225248154149067678929852793062368864372526397575528042259222936581665139 +/- 7.77e-77.

Thus density <0.510522525.

Updated baseline: arXiv:2512.04829v2 Table 1, https://arxiv.org/html/2512.04829, reports 0.5122645. Relative to its additive gap above pi^2/(15 sqrt(2)):
- Precise enclosure closes 3.705787188259372%.
- Outward-rounded 0.510522525 closes 3.705786795582668%.

The comparison uses the displayed newer baseline, not an undisplayed exact value. This does not attain the requested 20% gap closure.

Scope: source review, independent snapshot replay and objective reconstruction; not clean-room reimplementation, formal proof-assistant verification or independent reconstruction of external count certificates. No unrestricted literature-priority claim follows from this audit alone.

The original isolated snapshot is preserved as `audited_snapshot.json`. The portable `witness.json` changes only descriptive source/scope metadata; the runner asserts equality of all remaining fields. The independent objective enclosure above and the portable replay enclosure (radius below 2e-76) use different arithmetic orderings and both prove the same strict bound.
