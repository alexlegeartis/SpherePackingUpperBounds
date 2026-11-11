# Independent audit of the improved 48D cyclic certificate

PASS. The exact normalized objective is 18.148636519804413... <18.149. A safely rounded density consequence is Delta_48 <4.20664e-7. No optimality or literature-priority claim is established.

## Angular polynomial and scale lemma

Independent explicit Gegenbauer expansion exactly reconstructs the stored polynomial. The 24 coefficient slots end with four zeros; the actual degree is 19. Every coefficient is nonnegative, b0=1, and Q(1) matches B. Rational Sturm checks prove Q+10^-5<0 on [-1,1/2], Q<0 on [-1,5001/10000], and exactly one root between 5001/10000 and 1, with Q(1)>0.

For sorted triangles, z1<=z2 because the numerator difference is (u-1)(1+t)>=0. The single sign transition therefore excludes A=Q(z1)>0>C=Q(z2). For 0<=p<=q, F(T)=A exp(-pT)+C exp(-qT) satisfies F(T)<=max(F(1),0) for T>=1. Both same-sign cases are immediate. If A<=0<=C, the bracket C exp(-(q-p)T)+A is nonincreasing, and F has nonpositive derivative wherever it is positive. Thus it cannot cross from nonpositive to positive, and cannot exceed its initial positive value. This includes p=0 and p=q. The lemma correctly extends the triangle sign to all larger shortest-side scales and all damping parameters alpha>=64.

## Complete triangle verification

The independent verifier directly rebuilds exact Bernstein bounds instead of importing the producer's subdivision helpers. At 192-bit precision it certifies 3,238 rectangles and discards 463 rectangles solely because they violate the sorted-triangle domain. Maximum depth is 15. Signed secondary polynomial values are enclosed by centered polynomial expansions and rigorous exponential intervals. The compact domain is completely covered.

The analytic tails also pass: for t<=0 both exponent gaps are at least 1. For u>=5, a positive first secondary term forces gap>2*c*u-1>=4.001, while the other gap is at least 24. The old damping-36 text's 4.01 must not be copied with the new cutoff .5001; 4.001 is correct and sufficient. The exact inequality 64^30/30!>2B*10^5 proves the tail estimate. Collinear cases are included.

## Radial function and full bound

Independent reconstruction matches all 112 Fourier coefficients and verifies positive normalization and H(0)=1. Complete validated root isolation at 512-bit precision proves global Fourier positivity and physical negativity beyond squared radius 7. A separate 110-term Machin enclosure proves the required pi inequalities.

The new kappa=7/2 gives h=exp(-(7/2)*pi*(|x|^2-6)) and actual triangle damping 21*pi>64. The coupled radial condition after division by h is
P(s) exp((5/2)s-21*pi) <= M(2-exp(-(7/2)*(s-6*pi))).
Independent interval verification certifies all 2,105 intervals covering [6*pi,7*pi], with maximum depth 21. Beyond squared radius 7 the physical function is negative and the right side positive. Exact M and objective match the producer certificate.

The previously audited smooth PSD origin-feature construction is unchanged: origin value MB, two cross terms -Mh, and diagonal Mh^2. All hypotheses of the cyclic Cohn–de Laat–Salmon Theorem 2.4 with R=infinity are therefore satisfied. The full objective is P(0)+BM.

Using 18149/1000 and pi<355/113 gives the exact density upper bound (18149/1000)*((355/113)*(3/2))^24/24!, approximately 4.20663738872636e-7 and strictly below 4.20664e-7. This improves the previous explicit audited certificate but does not identify the optimal CE value or establish a literature record. The existing user-facing Page was not changed by this audit.

Artifacts: independent_cyclic64_angular.py/.json, independent_cyclic64_domain.py/.json, independent_cyclic64_radial.py/.json.
