# Annular mod-two coset cap (dimension n>=2)

Let Lambda have minimum norm >=1. In any coset of Lambda/2Lambda, distinct vectors differ by a vector of norm >=2; the coset is closed under negation. Consider representatives of distinct unoriented pairs in an annulus of squared norms [L,R], with R-L<=4.

For any two such representatives x,y, both |x-y| and |x+y| are >=2. Writing a=|x|²,b=|y|² gives |cos(x,y)| <= (a+b-4)/(2 sqrt(ab)). The last expression has partial derivatives with signs a-b+4 and b-a+4, so is nondecreasing in both variables throughout this annulus. Therefore |cos|<=s=1-2/R. For R<2 this is impossible for two pairs, so the cap is1.

For R>=2, choose integer k>=0 and Q(t)=t^(2k)(t²-s²). The spherical moment mu_{2k}=product_{j=0}^{k-1}(2j+1)/(n+2j) gives constant Gegenbauer coefficient Q0=mu_{2k}[(2k+1)/(n+2k)-s²]. Require Q0>0.

Every even monomial has nonnegative normalized Gegenbauer coefficients. For the coefficient of degree2j (0<=j<=k), the ratio between coefficients in t^(2k+2) and t^(2k) is
 (2k+2)(2k+1)/[(2k+2-2j)(n+2k+2j)].
The denominator decreases with j for n>=2, so this ratio is at least its j=0 value (2k+1)/(n+2k). Thus every coefficient of Q is nonnegative when Q0>0. Its new top-degree coefficient is positive as well.

For M representatives, positive definiteness gives sum_ij Q(cos_ij)>=Q0 M². Off-diagonal terms are nonpositive and diagonal terms equal Q(1)=1-s². Hence M<=Q(1)/Q0. Taking the floor and minimizing over any finite list of admissible k gives a rigorous cap, without requiring that the minimum over all k be found.

For a nonnegative radial majorant of height g_i on disjoint annuli, there are exactly2^n cosets. Thus its sum over all nonzero lattice vectors is <=2^(n+1) sum_i C_i g_i. This includes the zero coset and is deliberately conservative. Boundary points may be counted twice without harming an upper bound.

In the50D certificate, x=22|v|² and every original bin has width1/4 in x. Thus squared-norm width1/88<4. All cap computations are rational with k=0..500. The n=50 Gegenbauer positivity and prior rank-two sign theorem are separate dependencies.
