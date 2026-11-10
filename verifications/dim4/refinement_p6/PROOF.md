# Radial-degree-six four-dimensional refinement

Use the complete argument and four external neighbor-count assumptions in the [parent proof](../PROOF.md), with w(r)=(T_0(z),...,T_6(z)). Cutoff R=333/200, angular degree 16 and radial Fourier degree 31 are unchanged. Frozen finite-decimal rationals and exact matrix entries here define the functions; solver status does not enter the proof.

The unrestricted packing density is strictly below **0.633430327**. The reconstructed Arb enclosure is centered at 0.6334303260113454414754145742596239771574365762937475831059441425951904455187 with radius below 1e-76. This closes 13.881022% of the gap from the updated paper-reported baseline 0.63610277 to the D4 lattice lower pi²/16. It does not reach 20%, and comparison to a displayed eight-decimal baseline does not establish a worldwide-priority claim.

## Full-domain verification

All 17 matrices are exactly positive definite. Exact rational Sturm sequences prove global Fourier positivity and physical tail negativity. The count-adjusted pair inequality passes Arb subdivision in 215 boxes, depth 8. The far sign passes 17,723 boxes, depth 37; the cyclic sign passes 135 boxes, depth 13. All stacks exhaust. Independent snapshot replay reproduces these results; see the audit.

The far verifier now uses the exact cosine interval [-1,1/2]. This covers every required triangle with 1<=r,s<=R and opposite side d>=R: cos(theta)<=(r²+s²-R²)/(2rs)<=1/2. The last inequality follows because r²+s²-rs is separately convex, so its maximum on the radial rectangle occurs at a corner, where its values are 1, R², or 1+R²-R, all at most R². The endpoint is included. Thus the smaller cap changes computational cost, not the mathematical conclusion.

The cyclic numerator (abc)^16 times the cyclic sum is checked on the full side-length cube modulo permutation symmetry. Together with the far sign and zero support for two long sides, this proves every triangle condition. Half-open count cutoffs, closed support boundaries, Fourier normalization and the final pi²/32 volume factor are handled exactly as in the parent proof.

## External and computational assumptions

The count inputs 24,25,26,27 are the rigorous published average-degree bounds of de Laat–de Oliveira Filho–Vallentin, Table 2, https://arxiv.org/pdf/1206.2608, applied at slightly smaller radii. Their original SDP certificates are not reproduced here. Standard Fourier identities, Poisson summation, the positive-definite Gegenbauer addition theorem, exact rational arithmetic, Arb interval arithmetic and correctly-rounded IEEE-754/NumPy operations are trusted. This package is a reproducible computer-assisted certificate, not a proof-assistant formalization.
