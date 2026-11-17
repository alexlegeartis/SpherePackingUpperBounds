#83D parity-dependent extension: primary complete, final search audit pending

Use the verified80D Bachoc-Nebe parent and verified depth certificates101,108,203. geometry83.json selects three half-classes a,b,c. Full integral80×80 isometries derived from explicit coordinate permutations transfer a,b from101 and c from108. Two further exact isometries transfer a+c and b+c from203. build83.py verifies every full isometry, unimodularity, transported lift and residue, not merely a norm equality.

The remaining classes are301=a+b and302=a+b+c modulo2L. Both numerator norms are divisible by4, hence their half-coset distances are integers. The same integer-only code backend, with cutoff384 (distance3), completely excludes301 through3 in7,416,859nodes/0leaves. Cutoff512 (distance4) excludes302 through4 in149,419,024nodes/0leaves. Both runs complete all16glues. Thus their depths are at least4 and5, respectively, pending independent traversal. The smaller cutoffs retain the already audited local-shell completeness.

Vertical Gram is V=[[2,0,3/4],[0,2,3/4],[3/4,3/4,2]]. Its determinant is23/4. Exact minima by coefficient parity, indexed binary(a,b,c), are:
000(nonzero):8;100:2;010:2;110:4;001:2;101:5/2;011:5/2;111:3.
A finite exact check of[-3,3]^3 proves these globally: V−(9/10)I is positive definite, and any vector outside this cube has squared coordinate length≥16, hence vertical norm>8.

Horizontal depth bounds in the same order are0,6,6,4,6,11/2,11/2,5. Every nonzero parity sum therefore has total norm≥8, nonzero even vertical vectors already have norm≥8, and the embedded parent attains8. The Schur determinant is23/4 and exact center-density square is2^85/23. candidate83.json contains the full rational Gram, all isometries and all comparison fractions. A rational Machin enclosure gives ordinary density≥2.56962386473426e-18,16.5631010300% above the pinned preceding repository result. ExternalMW80×A3 is an attributable derived example, not asserted strongest.

The independent algebraic geometry checker verifies all maps and the vertical parity theorem using a different exact spectral bound and search cube. Final minimum certification still requires independent replay of301/302, not merely this geometric check.
