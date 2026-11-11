# 46D radial perturbation with inherited cyclic angular kernel

## Result and scope

The frozen exact certificate proves

**Delta_46 < 0.000000914882847.**

The complete primary verification gives the density enclosure centered at
9.1488284600262300906264110774671012774815093433564050608912408891658744426658146413084020013286703491110482354621588645003319793508420592842195065235930e-7,
with radius below 4.2e-159 at 512-bit precision. Relative to the pinned Cohn-table comparator 0.0000009295151556 this is about 1.5742% smaller. The comparator's live page returned 502 on 2026-11-11; the indexed table was crawled four months earlier. It cites Afkhami-Jeddi et al., arXiv:2006.02560v2. This is not an exhaustive worldwide-priority claim.

The source angular certificate and cyclic proof are those in alexlegeartis/SpherePackingUpperBounds, commit e756587e8fdf7827f085d6d0372cac267b110bd1, copied under `continuum_three_point/`. The radial function is new and must be verified in dimension 46; Fourier positivity is **not** inferred from the old 48D radial function.

## Exact witness

`certificate.json` defines rational polynomials P and H of degree 111 and a positive rational M. With X=pi|x|^2, set

f(x)=P(X)exp(-X),
h(x)=exp(-(17/8)pi(|x|^2-6)).

The physical minimum center separation is sqrt(6). The verifier reconstructs every Fourier coefficient from P using the dimension-46 Gaussian Fourier identity. It checks exact equality with H and H(0)>0. Complete certified root isolation establishes H>0 on [0,infinity), using no forbidden real roots and its positive value at zero. Thus the dimension-46 Fourier transform is nonnegative everywhere.

Complete certified root isolation and an exact endpoint sign establish P<0 on [7pi,infinity). On 6<=|x|^2<=7 the verifier proves

f(x) <= M(2h(x)-h(x)^2).

It uses rational alternating-series Machin bounds for pi, centered Arb polynomial evaluation, monotone exponential bounds, and complete interval subdivision. All 35 intervals pass, maximum depth 22. For |x|^2>=7, P<0 while 0<h<=1, so the coupled inequality follows there as well. These are global proof checks, not sample-point constraints.

## Angular kernel and optimized origin

Let Q be the exact degree-19 polynomial from `cyclic_angular40_exact.json`, and B=Q(1). Reexpand this same polynomial in normalized dimension-46 Gegenbauer polynomials:

Q(t)=sum_l b_l P_l^(46)(t).

Every coefficient is verified nonnegative over the rationals; sum_l b_l=B. Its positive constant coefficient is b0=1.8456112044750994..., independently checked from spherical moments. For both nonzero arguments at allowed separations, set

K(x,y)=(M/B)h(x)h(y)Q(xhat dot yhat).

Use degree-zero origin feature -sqrt(MB/b0), with nonorigin degree-zero feature sqrt(Mb0/B)h. Higher angular features vanish at the origin. Consequently

K(0,0)=MB/b0,
K(x,0)=K(0,x)=-Mh(x),
K(x,x)=Mh(x)^2.

The exact rank-one Gram identity (MB/b0)(Mb0/B)=M^2 proves the scalar origin extension is PSD; the nonnegative harmonic coefficients and addition theorem prove the full kernel is PSD. A smooth radial cutoff chi, zero near 0 and one for |x|>=sqrt(6), extends all features smoothly and rapidly decreasingly. In degree zero use sqrt(Mb0/B)chi*h-sqrt(MB/b0)(1-chi), and in positive degrees use chi*h times spherical harmonic features. This preserves every required value.

## Full cyclic domain

The angular polynomial, Gaussian damping, and squared minimum separation are unchanged from the source48 cyclic proof. The actual cyclic expression depends only on triangle side lengths, not on the ambient dimension. Multiplying by the new positive scalar M does not change its sign. The source's complete compact-domain proof and analytic tails therefore apply unchanged: alpha=51pi/4>40, 11,932 compact rectangles, 445 exact infeasibility exclusions, plus the scale and tail arguments. Both source angular and cyclic-domain verifiers are replayed by the entrypoint.

## Density conclusion and exact rounding

The existing Cohn–de Laat–Salmon cyclic continuum theorem, or the explicitly stated periodic PSD/cyclic summation proof in the source notes, applies. This is the cyclic generalization, not the literal individual-sign premise of Theorem 1.1. Normalizing f and K by H(0) gives

Delta_46 <= pi^23 (3/2)^23 / 23! * (P(0)+MB/b0)/H(0).

The verifier checks the rational objective, encloses the density with Arb, and separately proves strict rounding by substituting the rational upper Machin bound for pi. All sums converge absolutely because the features are rapidly decreasing. Periodic approximation extends the conclusion to unrestricted packings, not just lattices.

## Discovery and rejected candidates

Discovery began with t times the old normalized degree-111 polynomial plus a degree-15 linear combination of normalized dimension-46 Laguerre Fourier eigenfunctions. A floating-point LP imposed sampled sign and coupling conditions. Its first apparent improvement was invalid: the exact audit found ten forbidden Fourier roots and fourteen forbidden physical roots. Seven adaptive rounds added exact root-isolation counterexamples to the search, after which the final rationalized witness passed all global conditions. A one-part-per-million upward factor on M provides certified coupling slack; it is included in the exact stored M and final objective.

No optimizer status, grid check, or unsaved correction is used in the proof. The frozen exact coefficients are sufficient for verification. This computer-assisted proof is not a proof-assistant formalization or external peer review.
