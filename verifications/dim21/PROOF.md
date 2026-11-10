# Dimension 21 continuum density certificate

Producer and independent verification both pass. The unrestricted packing density satisfies Delta_21 < 0.005149554. The internal objective evaluates to density 0.00514955314370872549..., but the rounded outward decimal is used in comparisons.

Use Cohn–de Laat–Salmon, Three-point bounds for sphere packing, Theorem 1.1 (https://arxiv.org/abs/2206.15373), with r²=21/8 and R²=63/20. Five equally spaced radial hat nodes and angular degrees 0 through 16 define the same noncyclic kernel construction as the repository's 19D and 20D certificates.

The exact witness specifies rational P,H and matrices A_0,...,A_16. Set N=H(0)>0, f2(x)=P(pi|x|²) exp(-pi|x|²)/N. Its Fourier transform is H(pi|xi|²) exp(-pi|xi|²)/N. The radial Laguerre parameter is 19/2, with normalized basis L_k^(19/2)(2X)/binomial(k+19/2,k) and Fourier eigenvalues (-1)^k. Both direct monomial transformation and exact Laguerre reconstruction verify this identity.

Let w(u) be the five radial hat functions on [r²,R²]. Define phi_0(x)=(0,w(|x|²)) in the annulus, (1,0,...,0) at the origin, and zero elsewhere. The degree-zero term is phi_0(x)^T A_0 phi_0(y). In the annulus add sum_l=1^16 (w(|x|²)^T A_l w(|y|²)) P_l^(21)(xhat dot yhat), with P_l(1)=1. Divide by N. The kernel is zero outside this support. The theorem does not require kernel continuity. Exact LDL establishes positive definiteness of all17 matrices; the Gegenbauer addition theorem proves global kernel PSD.

Exact Sturm chains prove H positive on [0,infinity) and P negative beyond R²*pi. For every one of561 unordered fine radial node pairs, an exact angular Sturm check proves negativity over [-1,c], where c bounds the admissible cosine throughout adjacent radial rectangles. Each coarse cell has eight subdivisions. Cosine monotonicity follows from 2r²-R²=21/10>0. Bilinear radial interpolation extends all corner checks to the full continuous off-diagonal domain.

On each radial cell the coupled expression is P(pi*u) exp(-pi*u)+q_i(u), with q_i an exact quadratic polynomial. Producer Arb Taylor verification covers the entire annulus in96 boxes. Independent verification covers it in50 terminal intervals; both have maximum depth7. Thus all continuum conditions, not merely sampled inequalities, hold.

The exact objective is (P(0)+A_0[0,0])/N, stored in certificate.json. The volume multiplier is 2*11!/22!*pi^10*(21/8)^(21/2). The independent proof squares this positive expression and substitutes a rational Machin upper bound for pi; an exact rational comparison proves the strict stated decimal upper bound.

Run `python verify.py` and `python independent_verify.py`. Both fail by assertion on any failed condition. The numerical solver status is not used as proof. General-packing benchmark: https://cohn.mit.edu/sphere-packing/ gives upper0.0051596603948176 and lower0.0029298291067157... . The rounded certificate decreases the upper by about0.196% and removes about0.453% of that additive gap, not20%. No worldwide priority claim is made.
