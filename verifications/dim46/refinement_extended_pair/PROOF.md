# Extended pair interval in dimension 46

Let s=6, kappa=17/8 and h(x)=exp[-kappa*pi*(|x|²-s)]. The angular Q is the same exact dimension 46 polynomial from the prior border certificate, with constant normalized harmonic coefficient 1, nonnegative coefficients and B=Q(1). The complete angular Sturm, compact cyclic, scale-extension and analytic-tail proofs are replayed unchanged.

For nonzero allowed vectors, K(x,y)=(M/B)h(x)h(y)Q(xhat dot yhat). Its origin extension has K00=MB, K0x=-Mh(x), Kxx=Mh(x)² and an exact rank-one degree-zero Gram identity. The kernel is not truncated at radius sqrt7. A smooth cutoff is used only near zero; all features decay rapidly at infinity.

The cyclic continuum argument requires f(x)+2K(0,x)+K(x,x)<=0 for every |x|²>=6. Equivalently, with f(x)=P(X)e^-X and X=pi|x|², it requires f<=M(2h-h²). It does not require f to be nonpositive immediately after squared radius 7.

This certificate proves the coupled inequality on the complete interval X in[6pi,40], then proves P<0 on[40,infinity). Since 0<h<=1 for every allowed radius, the positive pair allowance M(2h-h²) makes the physical tail sufficient beyond 40. Complete Arb subdivision, using centered polynomial expansions and monotone exponentials, certifies 6,649 pair intervals with depth 28. No sampled-only pair point is used in the proof.

All 112 Fourier coefficients are reconstructed exactly in dimension 46. Complete certified root isolation proves H>0 on[0,infinity) and P<0 on[40,infinity), with exact endpoint signs. The degree 47 correction and exact100ppm guard in M are reconstructed independently by verify_membership.py.

The direct periodic PSD/cyclic summation theorem, followed by periodic approximation, therefore bounds unrestricted packing density by

pi^23 (3/2)^23 /23! * (P(0)+MB)/H(0).

A rational upper Machin bound for pi proves this is strictly below 0.0000008777229084. This uses the cyclic generalization of the Cohn-de Laat-Salmon continuum framework (arXiv:2206.15373), not the literal individual-sign premise of its Theorem 1.1.

The former degree 31/rate saturation certificates imposed physical negativity from 7pi and hence do not apply to this enlarged template. This is a deliberate escape from that scope, not a contradiction of the lower certificates.
