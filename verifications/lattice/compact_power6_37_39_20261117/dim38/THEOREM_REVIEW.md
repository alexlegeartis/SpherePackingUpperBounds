# Independent complete audit, 38D

**COMPLETE INDEPENDENT AUDIT PASS. Lattice-only density upper bound 1.23521549304e-5.**

Only frozen input SHA256 `2aaeab7c68312654ab2f8808661c7564703123ff5f63a95c540ca760abacaff4` is covered. The input snapshot in this directory is immutable and distinct from any subsequently refitted candidate.

Exact rational checks establish symmetry and positive definiteness of every final harmonic block, the origin Schur condition, exact Bernstein-to-monomial transformations, both Fourier-polynomial involution directions, a positive factorization of the Fourier polynomial, exact residual cap counts and the stated density normalization. The actual dimension is 38: the Laguerre initial polynomial is L1=38/2-x and recurrence is Lk=((2k+38/2-2-x)L(k-1)-(k+38/2-2)L(k-2))/k. The degree is 81 and contact scale is 17. The radial profile is 19/6-13x/102, cutoff 323/13, with powers 6,9,12.

Independent direct-cosine Sturm and quadratic Bernstein-control cap checks pass all 153 control polynomials. Fresh isolated scalar replay passes 653 annuli and 1799 interval nodes. Both tail proofs pass: all 81 physical roots are freshly enclosed below the real tail, and exact rational expansion of P(180+t) has all coefficients strictly negative. Recorded upper heights and exact mod-2 cap budget are checked again after interval replay. Even-dimensional normalization is rational; odd-dimensional normalization uses a rigorous Machin-formula lower bound for pi.

The geometric theorem is the compact harmonic rank-two plane decomposition reviewed in the completed 42/43/44 packages. Its indispensable geometric premise is nonpositivity of the normalized three-line sum on the full reduced-plane cube. Algebra and radial inequalities alone do not establish this premise. Arb/Jet interval arithmetic is shared with the reviewed producer implementation; exact linkage, cap algebra and translated-tail checks are separate.

Complete main geometry and independent coverage replay PASS: 217727 nodes / 108864 leaves at 384 bits, empty stack. Trace SHA256 `6639f31c2b2fca6b4cab9b5fdbeef984f526b09548fcc3468bf99dd270276b94`. This replay traversed the entire trace from the root without prefix resumption. The largest replayed leaf upper bound was approximately -0.0004505863632387657.
