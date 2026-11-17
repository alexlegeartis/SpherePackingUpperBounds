# Independent complete audit, 37D

**COMPLETE INDEPENDENT AUDIT PASS. Lattice-only density upper bound 1.8268708324e-5.**

Only frozen input SHA256 `78f99a960215ea1d0c44e2e93faf0cb1e06c4be1dba99cb14495433c2b09c550` is covered. The input snapshot in this directory is immutable and distinct from any subsequently refitted candidate.

Exact rational checks establish symmetry and positive definiteness of every final harmonic block, the origin Schur condition, exact Bernstein-to-monomial transformations, both Fourier-polynomial involution directions, a positive factorization of the Fourier polynomial, exact residual cap counts and the stated density normalization. The actual dimension is 37: the Laguerre initial polynomial is L1=37/2-x and recurrence is Lk=((2k+37/2-2-x)L(k-1)-(k+37/2-2)L(k-2))/k. The degree is 81 and contact scale is 17. The radial profile is 19/6-13x/102, cutoff 323/13, with powers 6,9,12.

Independent direct-cosine Sturm and quadratic Bernstein-control cap checks pass all 153 control polynomials. Fresh isolated scalar replay passes 653 annuli and 3381 interval nodes. Both tail proofs pass: all 81 physical roots are freshly enclosed below the real tail, and exact rational expansion of P(180+t) has all coefficients strictly negative. Recorded upper heights and exact mod-2 cap budget are checked again after interval replay. Even-dimensional normalization is rational; odd-dimensional normalization uses a rigorous Machin-formula lower bound for pi.

The geometric theorem is the compact harmonic rank-two plane decomposition reviewed in the completed 42/43/44 packages. Its indispensable geometric premise is nonpositivity of the normalized three-line sum on the full reduced-plane cube. Algebra and radial inequalities alone do not establish this premise. Arb/Jet interval arithmetic is shared with the reviewed producer implementation; exact linkage, cap algebra and translated-tail checks are separate.

Complete main geometry and independent coverage replay PASS: 216069 nodes / 108035 leaves at 384 bits, empty stack. Trace SHA256 `25a774562f44e2899695bb1101d856c1ee100fdf91f4079aa6470127ba349974`. The conservative resumed-prefix ceiling is not a measured worst-leaf margin.
