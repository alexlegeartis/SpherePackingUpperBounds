# Independent 22D tetrahedral auxiliary audit

PASS: at squared close-edge cutoff q=21/20, every finite unit-separated packing in dimension22 satisfies 3T <= B E, where E counts unordered close edges and T unordered all-close triangles, with

B = 208821822550866169/62500000000000 < 3341.14916082.

This is an auxiliary counting inequality, not a sphere-packing density bound.

## Independent computation

The exact polynomial is reconstructed by an explicit rising-factorial Gegenbauer expansion in transverse dimension21, independently of the source recurrence. All harmonic coefficients are nonnegative, c0=1, and B=Q(1). Exact rational Sturm variations and endpoint signs prove strict negativity on [-1,11/32]. The stored success flags are not used.

Two independently implemented interval replays cover the entire squared-edge cube [1,21/20]^6, including nonphysical edge tuples, without geometric exclusions. `verify_vertex.py` uses coordinates based at a vertex and completed21,735 leaves at maximum depth21, recorded in audit.log. `verify.py` uses equivalent midpoint coordinates for tighter interval dependency and completed3,090 leaves at depth17, recorded in midpoint_audit.log. Both use160-bit Arb, automatic interval derivatives and the multivariate mean-value bound, and both freshly verified all domains. Exact rational bisection gives complete closed coverage. Interval positivity of every transverse-radius denominator and finiteness are checked on each box. Floating-point scores select subdivision order only; no acceptance decision uses those scores.

The source sharp projection-cap proof was reviewed: at most one long edge contributes its opposite base-edge kernel, whose cosine is at most11/32; with two or more long edges every base-edge kernel vanishes. The strict far sign includes equality. The full compact cube covers all six close edges and every support boundary.

## Counting normalization

For an edge with N common neighbors, the positive kernel with origin BM, origin-neighbor -M and neighbor-neighbor (M/B)Q yields 0 <= BM-MN+(2M/B) sum_unordered_pairs Q. Summing over edges counts each all-close triangle exactly3 times and each tetrahedron once per supported base edge. The latter sum is nonpositive by the full-domain/far proofs, giving3T <= BE. There is no missing factor2 and no flooring of B. This is not an individual common-neighbor bound.

## Comparison and structural distinction

The prior sharp affine-envelope certificate was freshly replayed. Its exact B is1923216335019943159/500000000000000 <3846.43267004. The new exact constant is13.136418925559614% smaller; the prior certificate is not a proved optimum of its family.

There is nevertheless an exact separation from the affine family for this particular new polynomial. At t=3/8, Q(t)<=a(t-1/3) requires

a >= 7034984040475208405061/65536000000000000000 >107.34533753.

At t=17/59, it requires

a <= 1191611909750310730113443859/15728334578000000000000000 <75.76211607.

Both evaluations were independently reconstructed exactly. The incompatible inequalities prove that this Q cannot have any affine envelope on the sharp compact interval. Full joint tetrahedron constraints therefore strictly enlarge that feasible polynomial family. No claim about literature novelty, optimal auxiliary constants or density improvement follows.
