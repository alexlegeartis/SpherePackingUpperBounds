# Independent 5D early-tail spherical-lens audit

PASS: unrestricted sphere-packing density <0.509953498.

Frozen exact source SHA256:28585450c46b4fbdc6f9eefec55b7e6d29f18fc527b654d00359ed425cb1f949. Exact witness, prerequisites, source and embedded budget are copied here. This proof uses the conservative transverse projection formula, not the separately proposed sharper long-base cap.

## New geometric and spherical inputs

Write r=27/25, c for the long-base distance, and A=(r²-1)². Both incident distances of a common neighbor lie in[1,r]. Its axial coordinate relative to the midpoint has magnitude at most(r²-1)/(2c). Two neighbors therefore have axial difference squared at most D=A/c². Their transverse squared norms lie between L=1-c²/4-D/4 and U=r²-c²/4.

On c in[3/2,41/25], L is positive. The normalized cosine bound (u+v+D-1)/(2sqrt(uv)) is increasing in both squared norms because L-U+1-D>0. Thus the cosine is at most s(c)=1-(1-D)/(2U). Its derivative is negative: with x=c², the subtracted quotient has positive derivative numerator (x-A)²+A(4r²-A). All norm, monotonicity and denominator margins are checked rationally. Positive transverse norms and s<1 exclude zero vectors and duplicate normalized directions.

The three transverse code polynomials at c=3/2,31/20,8/5 are independently reconstructed in ambient dimension4 using normalized Chebyshev U_k/(k+1), rather than the source Gegenbauer recurrence. Nonnegative coefficients, constant coefficient1, exact global Sturm negativity on[-1,s], and exact Q(1) are checked. The addition theorem yields M²<=M Q(1), giving integer caps10,9,8. Their validity persists up to the next threshold by the proven angular monotonicity.

At c>=41/25, transverse pair dot products are bounded by U-(1-D)/2<0, decreasing with c. Strictly obtuse vectors in dimension4 have at most5 members: an affine dependence among6 vectors, split into positive and negative coefficients, would express a squared norm as a strictly negative cross sum. The zero-vector case permits at most one member and causes no exception. Taking the minimum with the independently recomputed centroid integer caps gives the later table5,5,5,3,2,1. No externally supplied spherical-code optimum is assumed; the specific polynomial certificates are fully replayed.

## Counting and endpoints

The previously independently reviewed lens incidence argument applies unchanged. Each eligible triangle has two sides<=27/25 and a unique long base>3/2. Both triangle and pair terms have multiplicity2 in the summed root-PSD/Fourier expansion. The verified f(d)+B(d)eta(d)<=0 funds the allowed positive far kernel with no new origin cost. Enlarging a periodic fundamental cell beyond the interaction range preserves all physical-center pair and triangle orbits; the incidence count and standard unrestricted-packing passage remain valid.

At exactly d=R=3/2, eta is zero, so the compact three-vertex condition remains unchanged. At every internal step endpoint eta is the smaller adjacent value; at the last support endpoint it is zero. The full far checker conservatively requires zero allowance for boxes meeting the exterior and checks both sides of the short-radius boundary. There is no double use of a bin budget.

## Independent numerical proof

Explicit half-integer Laguerre reconstruction/direct Fourier transform, all66 exact PSD pivots, global exact Sturm signs, Schur saturation, positive normalization and rational pi²/60 outward rounding pass. Both Gaussian correction weights are zero.

The density interval is centered at0.509953497415345158035896643175835425846568757645368535451747520236931752215657945308689001441897743162707534570006 with radius below2.86e-115.

Independent pair verification resolves1,042 closed leaves and explicitly checks all four published-count endpoints with the smaller allowance. The source-reviewed compact cyclic checker freshly passes99 boxes, including6 valid symmetry exclusions. It is not independently reimplemented here.

Independent budget verification uses derivatives with respect to X=pi*d² and proves all132 closed intervals using145 Taylor leaves, maximum depth4. It binds the physical coefficients, independently verified cap table and exact endpoint minima.

The independently written far checker uses explicit Gegenbauer polynomials, exact tensor Bernstein initialization, outward subdivision and independently derived exact distance extrema. It resolves2,879 boxes:1,011 passing leaves, of which113 use positive allowances, and430 valid geometry exclusions, maximum depth22. No uncovered box, numerical-grid inference or stored success flag is accepted.

## External assumptions and comparison

The original near-neighbor counts45,46,47,48 remain external published inputs from https://arxiv.org/html/1206.2608 Table2, with each threshold reduced by1e-6. Their original spherical SDP certificates are not rerun. The new lens polynomials, by contrast, are fully checked here.

The benchmark is the paper-reported .5122645 in https://arxiv.org/html/2512.04829v2 Table1; its witness is not replayed. D5 lower density is pi²/(15 sqrt2). This is a computer-assisted audit, not proof-assistant formalization or a worldwide-priority claim.
