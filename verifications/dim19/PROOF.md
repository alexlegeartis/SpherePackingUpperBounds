# Dimension 19 continuum density certificate

Producer and independent verification both pass. This is a general sphere-packing density bound, not a kissing-number or lattice-only bound.

Claim: Delta_19 < 0.009618351.
Exact objective: 278406330902078076451/9536743164059409037.
Density interval: 0.00961835093382645272416598000387883519167140758549870892621568166245471883098, with certified Arb radius below 4.35e-78.

Use Cohn–de Laat–Salmon, Three-point bounds for sphere packing, Theorem 1.1 (arXiv:2206.15373), exactly as in the independently audited 20D construction. Here r²=19/8, R²=57/20, with five equally spaced radial hat nodes. The kernel is zero outside the closed annulus except for its isolated origin coordinate. The theorem does not require kernel continuity.

The exact witness is certificate.json. It specifies P,H and rational matrices A_0,...,A_16. Set N=H(0)>0 and f2(x)=P(pi|x|²) exp(-pi|x|²)/N. Its Fourier transform is H(pi|xi|²) exp(-pi|xi|²)/N. The normalized radial Laguerre basis is L_k^(17/2)(2X)/binomial(k+17/2,k), with eigenvalues (-1)^k. The generalized binomial is computed exactly by successive multiplication by (19/2+k-1)/k.

For annular x, let w(|x|²) be the five radial hat functions, and phi0(x)=(0,w). At zero set phi0=(1,0,...,0), and elsewhere phi0=0. The unnormalized degree-zero kernel is phi0(x)^T A0 phi0(y). Add, on the annulus, sum over l=1,...,16 of (w(x)^T A_l w(y)) P_l^(19)(xhat dot yhat), where P_l(1)=1. Divide the total by N. All 17 matrices are positive definite by exact LDL; the addition theorem gives kernel PSD.

Exact Sturm chains prove H>0 on [0,infinity), and P<0 on [R²*pi,infinity), using a rational lower bound for pi. Every one of 561 unordered fine-grid corner angular polynomials is strictly negative on [-1,c], verified by exact Sturm chains. Each coarse radial cell is subdivided into eight pieces. The cosine upper bound c uses the upper neighbors of each corner and satisfies the defining squared inequality exactly. Monotonicity holds because 2r²-R²=19/10>0. Bilinear radial interpolation then proves the entire off-diagonal continuum sign, not just sampled points.

The coupled annular function is P(pi*u) exp(-pi*u)+q_i(u), with exact quadratic q_i on each coarse cell. The existing 256-bit Arb Taylor verifier proves this negative on the complete annulus, using 40 boxes and maximum depth 7.

The density multiplier is vol(B_(r/2)^19)=4^10*10!*pi^9/20! * (19/32)^(19/2). A 256-bit Arb enclosure gives the stated strict decimal upper bound.

Reproduction:
1. python verify.py
2. python independent_verify.py

The verification scripts print their complete pass summaries. The exact objective and rounded density bound are also recorded in certificate.json. Numerical precursor: snapshot19_final.json; its solver status is not used in the proof.

Benchmark: the general-packing table at https://cohn.mit.edu/sphere-packing/ lists 0.0098179551395438 and known lower density 0.0041473691458627... . The new upper value is about 2.03% smaller and closes about 3.52% of that gap. This is not a 20% gap improvement. No priority or novelty claim is made pending independent audit and literature review.
