# Optimized cyclic angular certificate at damping 36

Producer verification only until independently audited. This is a sign certificate, not by itself a density bound.

The polynomial in cyclic_angular36_exact.json is degree19, with nonnegative normalized dimension48 Gegenbauer coefficients and constant coefficient1. Exact Sturm tests establish Q<-10^-5 throughout[-1,1/2], Q<0 throughout[-1,501/1000], and exactly one root in(501/1000,1), with Q(1)>0. Thus Q changes sign exactly once on[-1,1]. Its value B=Q(1) is approximately969801987.975.

Use the sorted triangle variables t,u,v and cosines z1,z2 from CYCLIC96_PROOF.md. At shortest side a=1, the actual cyclic expression is

 Q(t)+A exp(-p)+C exp(-q),

where A=Q(z1), C=Q(z2), p=36(1-2ut), q=36(u²-2ut). Here0<=p<=q, and z1<=z2. Since Q has one sign change, A>0 implies C>=0.

For every T>=1, set F(T)=A exp(-pT)+C exp(-qT). If A,C>=0 then F decreases. If A,C<=0 then F<=0. The only remaining sign pattern is A<=0,C>=0. Wherever F(T)>0, its derivative is nonpositive, because C exp(-qT)>|A|exp(-pT) and q>=p. Consequently, in all cases,

 F(T)<=max(F(1),0).

Thus Q(t)<0 and Q(t)+F(1)<=0 imply Q(t)+F(T)<=0 for allT>=1. This includes all shortest-side lengths a>=1 and all Gaussian damping parameters alpha>=36 by T=a² alpha/36.

For0<=t<=1/2 and1<=u<=5 with2ut<=1, certify_cyclic_optimized_actual48.py verifies Q(t)+F(1)<0. It uses exact Bernstein upper bounds for Q(t), interval centered Taylor expansions for Q(z1),Q(z2), and Arb interval exponentials. Result cyclic_optimized_actual36_exact.json records11,199 successful rectangles, maximum dyadic levels19, precision128bits.

The remaining domains use the positive-part tail argument from CYCLIC96_PROOF.md. When t<=0 both exponent gaps are at least1. When u>=5 and t>=0, a positive first secondary term forces its gap>4.01, and the second gap is at least24. In either case, it suffices that 2B exp(-36)<10^-5. This is certified by the rational Taylor sum of exp(36) throughdegree60.

The angular data were obtained by a sampled LP, rationalized, and shifted downward by1/1000 before renormalizing the constant Gegenbauer coefficient. Only the final exact Sturm and interval checks, not the sampled optimization, establish the claims above. Kernel PSD and origin-feature extension are as in the earlier complete cyclic proof; a new radial function and coupling audit are still required for a new density claim.
