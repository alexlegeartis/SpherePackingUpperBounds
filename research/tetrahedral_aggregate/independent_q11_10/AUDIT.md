# Independent 22D q=11/10 tetrahedral auxiliary audit

PASS: 3T <= B E for finite unit-separated packings in dimension22, with unordered close edges defined by squared length <=11/10 and all-close triangles, where

B = 1271745209351204667/250000000000000 <5086.98083741.

The frozen certificate SHA256 is9e2bc0ad3fa7354fe0c0280309bb6d170fab8ee1f72a6197532a0759708aa0c9. This is an auxiliary inequality, not a density bound.

The independent implementation reconstructs normalized Gegenbauer polynomials in transverse dimension21 by explicit rising-factorial expansion, verifies c0=1, all harmonic coefficients nonnegative, B=Q(1), and exact outward rounding. A rational Sturm chain proves strict far negativity on [-1,6/17].

Full interval evaluation covers all six squared edges in [1,11/10]^6. The previously independently reviewed projection formulas, PSD origin normalization and3T counting factors are unchanged. Interval arithmetic checks positive transverse-radius denominators and finite gradients; exact rational subdivisions and closed faces preserve all boundaries.

The additional monotonicity rule is valid: a derivative of fixed sign on a closed box makes the maximizing endpoint dominate every point along that coordinate. The producer reduces every such coordinate jointly; this is sound because all sign bounds hold uniformly on the original box and therefore on each intermediate face. The independent implementation deliberately reduces only one coordinate at a time and recomputes all interval gradients. It freshly completes20,314 processed boxes,9,954 terminal leaves,407 endpoint-face reductions, maximum depth22, and zero geometric exclusions. Every face and box is resolved; no timeout or unverified remainder is accepted. Floating scores determine subdivision order only.

The prior sharp affine certificate was freshly replayed: B_old=160845191916055723/20000000000000 <8042.25959581. The new exact B is36.74687099058966% smaller. This compares two valid certificates, not proved optima.

The new polynomial also has an exact affine-family obstruction. At compact upper endpoint7/17, Q(t)<=a(t-1/3) requires

a >=24809255899410977889499233/265513259000000000000000 >93.43885873.

At lower endpoint7/29 it requires

a <=1375892903892437900868699/41022298000000000000000 <33.54012260.

Both rational evaluations were independently reconstructed. Their incompatibility proves this polynomial is outside the affine-envelope feasible family. No literature-novelty or packing-density claim follows.
