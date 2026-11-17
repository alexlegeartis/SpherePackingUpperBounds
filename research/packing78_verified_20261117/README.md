# Verified78D three-fiber Mordell–Weil packing

Ordinary packing density **>=1.41519803035180e-16**. This is a periodic three-translate construction, not a claimed78D lattice.

Method: primitive complex-multiplication dual normals of squared norm162/41, three compatible fiber projections of the MW80 parent, and exact coarea/shift geometry. The resulting minimum squared distance is3320/243, attained by an explicit source vector; the period Gram determinant is3^9*41^10. The source parent has minimum14 and determinant41^12.

Run `python -B run_all.py` with numpy and python-flint installed, without `-O`. The hash-checked wrapper uses a temporary copy. It replays the producer's parent construction, a separately implemented finite-field/height proof, the independent exact period/shift/density checks and the independent external comparison derivation. The separate parent implementation reconstructs190 polynomial sections and all17,955 pair heights using plain Python finite-field arithmetic; integral generation is proved by gcd1 of maximal minors. The external Shioda height theorem is explicitly required, as explained in packing78_parent_independent_20261117/AUDIT.md. Equality with the entire Mordell–Weil group is not assumed or claimed.

## Pinned external-derived comparison, not a strongest-prior assertion

[Bachoc–Nebe(1998), Extremal lattices of minimum8 related to the Mathieu groupM22](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf) and the [authors' L80 catalogue](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/L80.html) supply an even unimodular minimum8 parent. The first two catalogue basis vectors have Gram[[8,-4],[-4,8]]. They span a primitive normal2A2, so their perpendicular78D section has determinant48 and inherited minimum8. Its squared center density is2^78/48=2^74/3, ordinary density approximately9.523775961404e-17.

This is a pinned external-derived example, not a78D numerical record printed in the paper and not a proven strongest prior result. The new exact squared-density ratio gives an approximately48.5963168482% increase over this example. The separately pinned internal predecessor gives an approximately3.5823422851% increase; that internal value is not the principal comparison.

The catalogue Gram was independently parsed, its determinant and positive definiteness checked, and the comparison derived directly from its first two basis vectors. The source snapshot/hash, exact computation and result are included in packing78_independent_20261117. No exhaustive worldwide-priority claim is made.
