# Independent review: 50-dimensional lattice-only bound

Audited strict density upper: **1.75102693e-7**. This is a lattice-only bound, not a bound for all packings. No strongest-prior-source comparison or global-priority claim has been established.

## Mathematical review

The inherited scalar rank-two theorem is precisely the already published guarded degree-14 angular polynomial with profile h(t)=t^96 exp[-120(t-1)] and minimum length at least1. The exact coefficient source is compared with the pinned repository file at commit af863e2014561e060af7ee2b0eecffc468f4eca3. The independently coded exact triangular Gegenbauer re-expansion in dimension50 has all coefficients nonnegative and positive constant coefficient. Positive rescaling preserves the dimension-independent rank-two sign theorem. The old48D radial function is not used. This audit reuses that published, previously audited universal geometric proof rather than asserting a new replay of its entire compact gate.

For the new kernel, the constant-harmonic block is congruent to [[s,z],[z,a]], with a>0 and s=z²/a+10^-12. It is strictly PSD. Nonconstant harmonic blocks are nonnegative. On each primitive unoriented line, put H=sum_{k>=1} h(k²t). The successive ratios are at most exp(-168)<rho=10^-70, since 192log(1+1/k)-120t(2k+1)<=192/k-120(2k+1)<=-168. The exact factorial lower bound for exp(168) is checked. Hence H²<=(1+rho)/(1-rho) sum h_k². All even angular kernels have contact Q(1)=q. Summing independent pairs by saturated two-dimensional sublattices is legitimate by absolute convergence and gives a nonpositive contribution. Therefore positivity of the complete kernel sum implies 0<=s+4aH²/q+4zH/q summed over primitive lines. The corrected radial function is f+2aLh²/q+2zh/q, with L=(1+rho)/(1-rho).

For each mod-two coset and each squared-radius annulus [A,B] of width<=4, representatives of distinct unoriented pairs satisfy both |x-y|>=2 and |x+y|>=2. The function (a+b-4)/(2sqrt(ab)) has nonnegative partial derivatives on this square, giving |cos|<=1-2/B. If B<2, at most one pair is possible. Otherwise use t^(2k)(t²-c²), c=1-2/B. Its constant Gegenbauer coefficient is mu_2k[(2k+1)/(n+2k)-c²]. The ratio of the degree2j monomial coefficients is (2k+2)(2k+1)/[(2k+2-2j)(n+2k+2j)], minimized at j=0 for n>=2. Thus positive constant coefficient implies all coefficients nonnegative; off-diagonal values are nonpositive. The spherical PSD sum gives at most floor((1-c²)/constant) representatives. Finite optimization over k=0..500 is sufficient; global optimality is not required. Same-direction duplicates cannot occur in these annuli under the same separation inequalities.

Exactly2^50 cosets, including zero, and both signs yield the majorant budget 2^51 sum cap_i*height_i. Counting shared annular boundaries twice is conservative. All632 original bins cover x=22 through180 without gaps. Beyond180, P is negative by complete root isolation, and aL+z<0 together with0<h<=1 makes the kernel correction nonpositive. Combining the radial majorant with PSD gives sum_{v!=0}f(v)<=s+budget.

With x=pi|v|² and contact x=22, the independently checked Fourier normalization is Q_F(0)=1. Scaling to minimum length1 makes the integral (pi/22)^25. Poisson summation therefore gives density <=(22/4)^25/25! [P(0)+s+budget]. The final guard is10^-12, NOT the earlier preliminary normalization script's10^-20. The independently reconstructed rational expression is strictly below1.75102693e-7.

## Computational review

- Exact factor reconstruction verifies all18 rational double-root factors,17 positive quadratics and one positive linear factor on the nonnegative ray. Unknown factor types are rejected.
- A separately coded direct Gaussian monomial transform matches all72 Fourier coefficients; it does not rely on the producer's Laguerre recurrence.
- Complete1024-bit root isolation counts all71 spatial roots, with33 possibly real roots below180 and P(180)<0. Thus no spatial positive window beyond the cutoff is missed.
- Reviewed combined degree212 Taylor polynomial and exponential-remainder bound: each shifted polynomial multiplies the degree20 exponential Taylor polynomial; the remainder is bounded by its full-bin absolute maximum times exp(lambda*r)(lambda*r)^21/21!. Outward interval enclosures include every term.
- Replayed the final annular implementation at1024bits, adding exact adjacency/domain checks and binary-forest identity nodes=2*leaves-632. All632 bins pass,1612 total nodes. Independent exact script recomputes every rational cap and budget, Schur guard, line constant and final rounding.

This is independent mathematical/code review and higher-precision replay of the annular implementation plus a separately written algebra/normalization checker. It shares FLINT/Arb with the producer; it is not a formal proof or an independent interval arithmetic backend. Original producer certificate is preserved separately; regenerated interval heights may differ in low outward-rounded bits.

## External comparison scope

Primary Henry Cohn/MIT sphere-packing density-table metadata lists dimensions1–48,56,64,72, not50: https://dspace.mit.edu/entities/product/d8a88a76-3fa7-4b75-89ef-ea71a0e05b77 and https://cohn.mit.edu/sphere-packing/ . A sampled CE objective is not a published or certified baseline. Afkhami-Jeddi et al., https://arxiv.org/html/2006.02560v2 , describes numerical bootstrap spectra; as noted in its Section2.3, lifting an optimal minus-one eigenfunction to a same-radius Cohn–Elkies auxiliary function is conjectural in general, and positivity certification is separate. No percentage improvement or SOTA assertion is justified by this audit alone.
