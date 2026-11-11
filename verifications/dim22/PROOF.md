# A Gaussian cyclic density certificate in dimension 22

## Claim

The exact witnesses and complete checks in this directory prove

**Delta_22 < 0.003653383161.**

The density enclosure is [0.0036533831609501341267126396728558645797789028984200797508912659434683062291395926515097012982093770455911300489588974235826541277608986328834876355392171 +/- 1.35e-155].

The pinned comparator is 0.0037259419689206 from Henry Cohn's sphere-packing table, associated with Afkhami-Jeddi et al., arXiv:2006.02560v2. The indexed table was retrieved on 2026-11-11 and had been crawled four months earlier; the live page returned 502. This is a specific published comparison, not an exhaustive worldwide-priority claim. The newer low-dimensional table in arXiv:2512.04829v2 must be used in dimensions 4–7 and 9–16; it is not substituted for this 22D comparison.

This is a computer-assisted upper bound for arbitrary packings, not a construction, kissing-number bound, optimality proof, external peer review or proof-assistant formalization.

## Exact functions and positive semidefinite kernel

Set s=11/4, kappa=51/11, X=pi|x|² and alpha=kappa*pi*s=51*pi/4>40. The rational degree 111 polynomials P,H in certificate.json satisfy

 f(x)=P(X)exp(-X),  fhat(y)=H(pi|y|²)exp(-pi|y|²)

with the Fourier convention exp(-2*pi*i*x·y). The verifier independently reconstructs all 112 coefficients of H in physical dimension 22. Complete certified root isolation covers all 111 roots with multiplicities, proves H>0 on[0,infinity), and proves P<0 on[7*pi,infinity). H(0)>0 is exact.

The actual degree 14 polynomial Q in continuum_three_point/cyclic_angular40_exact.json is expanded in normalized dimension 22 Gegenbauer polynomials P_l^n with P_l^n(1)=1. Its constant coefficient b0 is exactly 1 and every harmonic coefficient is nonnegative. Its value B=Q(1) is approximately 86892.28000427491; the stored rational is authoritative.

For nonzero vectors whose squared radii are at least s, define

 h(r)=exp[-kappa*pi*(r²-s)],
 K(x,y)=(M/B)h(|x|)h(|y|)Q(xhat·yhat),
 K(x,0)=-M h(|x|), K(0,0)=MB,

where M>0 is the exact rational in certificate.json. The normalized harmonic addition theorem gives PSD angular features. In degree zero, the origin feature is -sqrt(MB) and the nonorigin feature is sqrt(M/B)h; hence the origin Gram determinant is exactly zero and (MB)(M/B)=M². Higher-degree origin features vanish.

Choose a smooth radial cutoff chi that vanishes near 0 and equals 1 for r≥sqrt(s). Use sqrt(M/B)chi*h-sqrt(MB)(1-chi) in degree zero, and chi*h times the positive-degree angular features. This defines globally smooth rapidly decreasing features with the same required values at 0 and outside the excluded ball. Thus K is a smooth, symmetric, PSD, rapidly decreasing kernel. No radial discontinuity or compact-support assumption is needed.

## Pair inequality and full radial signs

On s≤r²≤7, complete Arb subdivision proves

 P(pi*r²)exp(-pi*r²) ≤ M(2h(r)-h(r)²).

There are 255 accepted intervals, maximum depth 15. Each interval is covered using a centered polynomial expansion and monotone exponential bounds, with rational Machin enclosures for pi. This is not sampling. Since K(x,x)=M h(r)², this is exactly

 f(x)+2K(x,0)+K(x,x)≤0.

For r²≥7, P<0 and 0<h≤1 imply the same inequality. The number 7 is a sign-check threshold, not the support of K: the kernel remains nonzero and Gaussian-decaying at all large radii.

## Complete cyclic triangle sign

The exact angular verifier reconstructs Q independently, verifies its nonnegative harmonic coefficients, and uses rational Sturm sequences to prove

 - Q(t)<-1/100000 for -1≤t≤1/2;
 - Q(t)<0 for -1≤t≤5001/10000;
 - Q has exactly one root between5001/10000 and1, with Q(1)>0.

For a triangle, scale the shortest side to1, let u≥1 be the second side, t the included cosine between them, and w=sqrt(1+u²-2tu) the longest side. Therefore t≤1/2 and 2tu≤1. Put z1=(1-ut)/w and z2=(u-t)/w. After dividing the cyclic sum by its positive largest Gaussian factor, its sign is the sign of

 Q(t)+exp[-A(1-2tu)]Q(z1)+exp[-A(u²-2tu)]Q(z2),

where A=kappa*pi*(actual shortest side)²≥alpha>40.

The degree 14 polynomial is padded to degree 19 in the Bernstein bounds; degree elevation preserves enclosure. The full compact domain 0≤t≤1/2,1≤u≤5,2tu≤1 is covered by rational dyadic rectangles. At A=40 the Arb verifier accepts 2770 rectangles and excludes 319 rectangles solely by the exact contradiction 2tu>1; maximum depth 15. Interval images enclose both secondary cosines and exponentials, and exact Bernstein coefficients bound the first polynomial term. Every accepted rectangle proves strict negativity. Degenerate triangle boundaries are covered.

To extend from A=40 to larger A, set p=1-2tu and q=u²-2tu, so0≤p≤q, and note z1≤z2. The single-sign-change property of Q allows only the following sign patterns for (Q(z1),Q(z2)): both nonnegative, both nonpositive, or first nonpositive and second nonnegative. For either same-sign pattern, the secondary sum at larger A is bounded by its value at 40 or by 0. In the mixed case, the function a*exp(-p*A)+b*exp(-q*A), a≤0≤b, has no interior maximum: its derivative changes sign at most once, from negative to positive. Its supremum on[40,infinity) is therefore at the left endpoint or its nonpositive limiting value. Since Q(t)<0, compact negativity persists for all A≥40.

The omitted domains follow analytically. For t≤0, both exponential gaps are at least 1, and the addition theorem bounds |Q|≤B on[-1,1]. For u≥5 and t≥0, w≥u and z1≤1/u≤1/5, so Q(z1)<0; the second exponential gap is at least u²-1≥24. The exact check 40^30/30!>2B*100000 implies exp(40)>2B*100000, so the possible positive secondary terms are smaller than the cap margin 1/100000. Hence the cyclic sum is nonpositive for every packing triangle, including all unbounded radii. No global-tail sampling is used.

## Direct periodic summation and normalization

This uses the cyclic refinement of the Cohn–de Laat–Salmon continuum framework (arXiv:2206.15373), not the literal individual-sign premise of printed Theorem 1.1. Here is the complete summation argument.

Let a packing be periodic with N centers per fundamental cell of a lattice Lambda of volume V and squared minimum separation s. For every center p, PSD gives

 0≤sum over q,r in the whole packing of K(q-p,r-p).

The inequality follows first for finite sums and then by absolute convergence. Sum over the N representative centers p. Terms with q=r=p contribute N*K(0,0). Terms with precisely two coincident centers contribute the ordered-pair sum of D(x)=2K(x,0)+K(x,x). Terms with three distinct centers can be grouped, modulo lattice translations, into twice the cyclic sum for each unordered triangle. Every such contribution is nonpositive by the previous section. Thus

 sum over ordered nonidentical pairs per cell of D(q-p) ≥ -N*K(0,0).

The pair inequality f+D≤0 gives the upper estimate

 sum over all ordered pairs per cell of f(q-p) ≤ N*(f(0)+K(0,0)).

Poisson summation and fhat≥0 give the lower estimate N²*fhat(0)/V for that same sum. Hence the center density obeys N/V≤(f(0)+K(0,0))/fhat(0). Multiplying by the volume of one packing ball of radius sqrt(s)/2 yields

 Delta_22 ≤ pi^(22/2)*(s/4)^(22/2)/Gamma(22/2+1) * (P(0)+MB)/H(0).

All sums converge absolutely because the features and f decrease rapidly. Grouping triangles therefore introduces no limiting sign assumption. Arbitrary packing densities are approximated by finite-window periodizations with a boundary gap, so the bound is not restricted to lattice or periodic packings.

The verifier checks the exact rational objective (P(0)+MB)/H(0), then uses a rational Machin upper bound for pi and, in odd dimensions, an exact rational upper bound for sqrt(s). This proves the final strict decimal independently of the printed Arb interval.

## Discovery and reproducibility

This is the first fresh optimization of this noncompact Gaussian-angular family in dimension 22 in this research campaign. It differs from the previously unsuccessful compact solid-hat runs. The angular LP uses fresh dimension-specific normalized harmonics. Rational export subtracts 1/1000 from the constant harmonic coefficient and renormalizes b0 to1; only the complete subsequent checks prove feasibility.

The physical polynomial was discovered as a scaled degree 111 source plus degree 31 Laguerre corrections in dimension 22. Initial LP outputs had forbidden Fourier and physical roots. Six exact-root-guided repair rounds produced this frozen witness. M has an explicit upward factor 501/500 for pair certification. These are discovery choices, not assumptions in the proof.

The standalone replay reads only the frozen rational radial and angular coefficients. It does not require the old 44D source, optimizer, search grids or any claimed near-feasible result. Run python verify.py with python-flint 0.9.0 and assertions enabled. Independent audit status is reported separately; filenames inherited from earlier independent verifiers do not themselves certify independent review of this new witness.
