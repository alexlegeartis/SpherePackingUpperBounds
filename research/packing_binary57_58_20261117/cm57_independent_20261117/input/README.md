# Frozen MW56 ramified coset, new57D/58D constructions

Status: primary full exclusion completed; independent audit pending. No publication authorization is implied.

## Parent and exact setup

Parent curve y²=x³+1+t30 overF29² with field generator²=2. The140 saved unitary sections have x=−(at+b)^10,y=(ct+d)^15. parent/verify_parent.py replays every curve equation, pairwise polynomial-gcd height, integral generation, positive definiteness and determinant29^8. All30 finite singular fibers are irreducible typeII and infinity is smooth. The elliptic-surface height formula withχ=5 gives every nonzero section norm>=10; saved polynomial sections attain10. Published parent input: Schütt–Shioda, Mordell-Weil Lattices,2017 author draft Theorem13.35 p357, sections p359; https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf . The minimum theorem and height formula are external mathematical dependencies.

The source was independently rebuilt here from140 sections, not assumed to share the earlier repository's120-section coordinate basis. cvp_basis_actions.json gives exact unimodularT and reducedGramQ=TGTᵀ. The target comes from the ramified CM map M=3(I−CM)⁻¹; parent/primal_cm.json stores its exact identities. target_class.json is random class37 in that older reduced basis; old_to_new transports it to the enumeration basis. The actual proof only needs the saved integer targetr and its congruence properties, not random generation or symmetry completeness.

## Complete exclusion

random56_37_bits.txt pins56, cutoff54, targetr, integerQ, and binary64 bit-pattern enclosures of the exact LDL decomposition in random56_37_exact.json. The reported complete run is

 EXHAUSTIVE_NO_VECTOR; nodes3791478090; seconds457.378652744.

The cap was5000000000, not reached. The exact input SHA256 is0f892d0ce87ab91a1d0faf4f9b08a16c072ea9af748bf86baad61d2319d8af99; verifier source SHA256 is60ddae41065e477d11347b6b8893b43abb7dd8c075508d1523e66f3407360392.

Replay: g++ -O3 -std=c++17 -ffp-contract=off -fno-fast-math -frounding-math verify.cpp -o verify; ./verify random56_37_bits.txt 5000000000 fast bits . The full run takes about8minutes on the current computer. A fresh same-code replay is useful but is not an independent algorithm.

For every integerz, the numerator norm (r−3z)Q(r−3z)ᵀ is0mod6 because Q is even andrQrᵀ=0mod6. Therefore excluding all numerator norms<=54 proves the half-open gap up to60. random37_contact.json gives an exact norm60 vector in the class, so the coset squared minimum is exactly60/9=20/3. It also gives original-parent coordinates for compact extension Grams.

## Enumeration completeness and arithmetic assumptions

PROOF.md describes the original interval engine. The fast mode encloses the computed center with an absolute2^−10 margin. Inputs are restricted to dimension<=110, target coordinates0,1,2; lower-triangular coefficients within±16 with enclosure width<=2^−47; every integer coordinate is restricted to±2^20, with overflow/range violations aborting rather than pruning. The center-error estimate under these bounds is below2^−13, leaving a large safety margin. Child contributions and residual radii are rounded outward using nextafter. Every enclosing integer is visited; no pruning coefficients are used. Exact signed128-bit leaf norms are guarded by Gram-entry<=2^30 and residual-entry<=2^22, bounding accumulated magnitude below2^88. Nodes cap/nonfinite/range errors return NOT_CERTIFIED, never an empty-ball result.

Dependencies: IEEE754 binary64 round-to-nearest, correctly rounded basic operations/square root on the replay platform, and the stated non-fast-math compiler flags. Runtime checks rounding mode and IEC559. The independent audit must verify these inequalities and the exact LDL/bit enclosures, not simply trust the completed-run label. The40 small exact-brute-force controls were tested previously; they are not a replacement for auditing the56D input.

## New lattices

Write x for the saved contact divided by3. In57D adjoin(x,sqrt(10/3)). Classes±x have minimum20/3; nonzero multiples of3 in the layer have vertical norm>=30. Minimum is exactly10, attained in the parent. Determinant=(10/3)29^8.

In58D take triangular vertical basisa,b with Gram[[10/3,5/3],[5/3,10/3]] and adjoin(x,a),(-x,b). Ifu−v is nonzero mod3, the horizontal contribution is>=20/3 and vertical contribution>=10/3. Ifu−v=0mod3, thenu²+uv+v² is a positive multiple of3 for every nonzero(u,v), so its vertical norm is>=10. Minimum is exactly10. Determinant=(25/3)29^8. The full rational Grams are in extension_candidates.json.

Center-density squares are(5/2)^n divided by these determinants. Ordinary densities are approximately1.533071110212405e−11 in57D and5.024193656484039e−12 in58D; rigorous rounded floors are to be verified in the audit.

Relative to the previous repository constructions, the exact ratios are sqrt(21/20) and21/20. External comparisons must not use those preceding records:57D uses the derived MW56×A1 product, giving ratio sqrt3;58D uses the previously source-audited Γ72 rank14 Leech-block section with center-density square2^36/3. The stronger baseline review remains a separate priority question, not settled by this construction proof.
