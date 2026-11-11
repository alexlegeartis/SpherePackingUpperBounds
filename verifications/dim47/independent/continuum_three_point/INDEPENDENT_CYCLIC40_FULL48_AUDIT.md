# Independent audit of the 48D cyclic bound below 17.835

PASS. The exact normalized objective is 17.834291498319086... <17.835. The safely rounded rigorous density bound is Delta_48 <4.13385739314e-7.

## Angular polynomial and all triangles

The independent explicit Gegenbauer expansion reconstructs every coefficient of the degree-19 angular polynomial in cyclic_angular40_exact.json. All harmonic coefficients are nonnegative, b0=1, and Q(1)=B=5499674981494692543917/6060000000000. Trailing zero coefficient slots do not raise the degree.

Exact rational Sturm checks prove Q+10^-5<0 on [-1,1/2], Q<0 on [-1,5001/10000], and exactly one root between that cutoff and 1, with Q(1)>0. Therefore the single-sign-transition scale lemma independently checked in INDEPENDENT_CYCLIC64_FULL48_AUDIT.md applies without modification. The ordering z1<=z2 excludes the only problematic secondary-term sign pattern. For every scale T>=1, the two-exponential contribution is at most its positive part at T=1. No unjustified monotonicity of the entire signed sum is assumed.

The independent 192-bit verifier directly rebuilds exact Bernstein bounds and certifies all 11,932 compact-domain rectangles. It discards 445 rectangles only for violating the sorted-triangle restriction; maximum depth is 19. The signed secondary terms are bounded with centered polynomial interval expansions and rigorous exponentials. Separate tails cover t<=0 and u>=5. The exact inequality 40^30/30!>2B*10^5 suffices; the correct first-gap bound in the latter tail is 4.001. All separated triangles and their scale/damping extensions are covered, including collinear cases.

## Radial function and complete continuum theorem

The independent radial verifier reconstructs all 112 Fourier coefficients directly from the exact rational physical polynomial. Normalization is positive and H(0)=1. Complete validated root isolation at 512-bit precision proves the global Fourier sign and the physical tail sign beyond squared radius 7. Exact endpoint signs fix the sign on each ray. Independent 110-term Machin bounds enclose pi.

The new kappa is exactly 17/8. Its neighbor feature is h=exp(-(17/8)*pi*(|x|^2-6)); rescaling minimum distance sqrt(6) gives alpha=51*pi/4>40. The correctly normalized coupled inequality after division by h is
P(s) exp((9/8)s-(51/4)*pi) <= M*(2-exp(-(17/8)*(s-6*pi))).
Independent rigorous interval subdivision certifies all 33 intervals covering [6*pi,7*pi], maximum depth 20. Beyond squared radius 7, physical negativity and 0<h<=1 finish the unbounded constraint. Exact M and the complete objective P(0)+BM agree with the producer's certificate.

The globally smooth PSD feature construction, origin value MB, cross terms -Mh and annular diagonal Mh^2 are unchanged from the previously audited cyclic theorem application. All requirements of Cohn–de Laat–Salmon equation (2.4) and Theorem 2.4 with R=infinity are therefore satisfied. The result is a density upper bound, not merely a cyclic sign lemma or sampled numerical solution.

## Density, comparison, and presentation corrections

Using 17835/1000 and pi<355/113, the exact rational density bound is
(17835/1000)*((355/113)*(3/2))^24/24!.
It is strictly below 4.13385739314e-7; this decimal direction was checked exactly. Compared with the stated numerical table normalization 19.5401246309 and extremal benchmark 1, the additive gap reduction is about 9.20%, not 20%. This does not identify the optimal CE value or establish literature priority.

The producer proof names certify_cyclic_optimized_actual64.py in one place; the actual candidate-specific script is certify_cyclic_optimized_actual40.py. This is a documentation typo, not a mathematical failure. The independent verifier reads cyclic_angular40_exact.json directly. Existing user-facing Pages were not changed by this audit.

Artifacts: independent_cyclic40_angular.py/.json, independent_cyclic40_domain.py/.json, independent_cyclic40_radial.py/.json.
