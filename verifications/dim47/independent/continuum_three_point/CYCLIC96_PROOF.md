# Lower damping cyclic certificate in dimension 48

This is an improved global cyclic-sign certificate, not yet a density bound. Use the already audited degree-19 Q with cutoff c=501/1000 and nonnegative Gegenbauer coefficients summing to B. The earlier exact Sturm certificate gives Q(t)<-10^-12 for -1<=t<=1/2.

Sort a triangle's side lengths as a<=b<=d, with a>=1. Set t=(a²+b²-d²)/(2ab), u=b/a, v=d/a. Then -1<=t<=1/2, u>=1, 2ut<=1, and v=sqrt(1+u²-2ut). The other two cosines are z1=(1-ut)/v and z2=(u-t)/v. Divide the cyclic Gaussian sum by its first weight. It is bounded above by

 Q(t) + max(Q(z1),0) exp(-alpha(1-2ut))
      + max(Q(z2),0) exp(-alpha(u²-2ut)).

Indeed the actual exponent gaps contain an extra factor a²>=1; negative terms may be dropped and positive terms decrease with that factor. This upper bound is also nonincreasing with alpha.

For 0<=t<=1/2 and 1<=u<=5, certify_cyclic_triangle48.py with argument 96 verifies the displayed upper bound is negative. It uses exact rational Bernstein coefficients for Q on dyadic t intervals, Arb interval evaluation of the two other cosines and Q, and interval exponential upper bounds. Rectangles with 2tu>1 everywhere are discarded. The run certifies 6,712 rectangles, reaching dyadic depths 16 in each variable, at 128-bit precision. The angular polynomial is known nonpositive through c, so those positive parts can be set to zero. Otherwise the bound min(B,Q(interval).upper()) is used.

For t<=0, both exponent gaps are at least 1, so the bound is at most -10^-12+2B exp(-96)<0. This last inequality follows already from 96^30/30!>2B*10^12.

For u>=5 and t>=0, 2ut<=1 forces t<=1/10. If Q(z1)>0 then z1>c, giving 1-2ut=2vz1-1>2cv-1>=2cu-1>=4.01. The second gap u²-2ut>=24. Thus the same exponential estimate proves strict negativity. All triangles, including collinear ones, are covered.

Consequently the Gaussian cyclic sign holds for alpha=96 and every alpha>=96. The smooth-cutoff positive definite extension from CYCLIC_GAUSSIAN_LEMMA.md applies. The radial coupling and objective still need separate verification.

The earlier two-variable relaxed bound without the actual triangle restrictions failed; it is not used here. The provisional alpha=64 interval attempt also did not certify and is not claimed feasible.
