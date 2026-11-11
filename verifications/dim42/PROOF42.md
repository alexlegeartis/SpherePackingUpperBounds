# Dimension 42 with a dimension-specific cyclic angular polynomial

## Result

The frozen exact witness proves

**Delta_42 < 0.000003677844262.**

The complete primary proof gives a density enclosure centered at 3.6778442617374951152832281446353112713623079583389e-6, with radius below 6.1e-158 at 512-bit precision. It improves the pinned Cohn-table comparator 0.0000038780970907 by approximately 5.1636880%. The indexed table was crawled four months before the 2026-11-11 retrieval; its live page returned 502. It cites Afkhami-Jeddi et al., arXiv:2006.02560v2. This is a pinned comparison, not an exhaustive current-worldwide-priority claim.

The underlying density theorem is the existing Cohn–de Laat–Salmon cyclic continuum framework, arXiv:2206.15373. This is a computer-assisted density upper bound, not a packing construction, optimality proof, external peer review or proof-assistant formalization.

## What changes from the source48 model

Both the angular polynomial and the radial Fourier pair are new dimension-42 data. The old 48D sign checks are not substituted for their verification. Only the structural PSD construction and analytic cyclic scale/tail lemmas are inherited. The scale and tail arguments are stated below; no historical numerical witness is a proof input.

- `continuum_three_point/cyclic_angular40_exact.json` defines the new rational degree-18 polynomial Q, expanded in normalized dimension-42 Gegenbauer polynomials. Its constant coefficient is exactly 1; every coefficient is nonnegative; B=Q(1)=119501506.83202717... .
- `certificate.json` defines the new rational degree-111 polynomials P and H, positive rational M, and exact objective.
- Minimum squared center separation is 11/2; the physical tail threshold remains 7; kappa=51/22 and alpha=(11/2)*pi*kappa=51*pi/4>40. The minimum radius differs from the previous 45D witness; the normalization and coupled interval were changed accordingly.

## Exact angular and complete cyclic proof

The angular verifier independently constructs the normalized dimension-42 Gegenbauer basis from its explicit rising-factorial formula, checks the polynomial coefficients, verifies every coefficient is nonnegative and that their sum is B, and proves by exact rational Sturm sequences:

1. Q(t)<-1/100000 for every -1<=t<=1/2.
2. Q(t)<0 for every -1<=t<=5001/10000.
3. Q has exactly one root in (5001/10000,1), with Q(1)>0.

For sorted normalized triangle variables t,u satisfying 0<=t<=1/2, 1<=u<=5 and 2tu<=1, set

w=sqrt(1+u^2-2tu), z1=(1-ut)/w, z2=(u-t)/w.

The exact interval domain verifier proves

Q(t)+exp(-40(1-2ut))Q(z1)+exp(-40(u^2-2ut))Q(z2)<0

over the whole domain. It passes 5,313 rectangles, with 411 rectangles excluded by the exact incompatible-domain condition 2tu>1, maximum depth 15 and 192-bit Arb precision. Polynomial Bernstein bounds and centered interval expansions, not sampled values, certify the accepted rectangles.

The same single-sign-change scale lemma applies to this new Q. The secondary exponents satisfy 0<=p<=q and secondary angles z1<=z2. Their signs are either both nonnegative, both nonpositive, or first nonpositive/second nonnegative. In each case A exp(-pT)+C exp(-qT)<=max(A exp(-p)+C exp(-q),0) for T>=1. Thus the checked inequality and Q(t)<0 extend to every larger shortest side and every alpha>=40.

The remaining domains use the same analytic tail argument. For t<=0, both exponential gaps are at least 1. For u>=5, a positive first secondary term forces its gap greater than 4.001, and the second gap is at least 24. Nonnegative normalized harmonic coefficients imply |Q|<=B on [-1,1]. The exact rational check 40^30/30!>2B*100000 therefore makes these positive contributions smaller than the negative cap margin. These arguments depend on the newly verified sign properties and B, not on dimension 48 or its old polynomial.

## New dimension-42 radial proof

Let X=pi|x|^2 and

f(x)=P(X)exp(-X),
h(x)=exp(-(51/22)pi(|x|^2-11/2)).

The verifier independently reconstructs all 112 Fourier coefficients in dimension 42 and compares them exactly with H. Complete certified root isolation and an exact positive endpoint establish H>0 on the entire nonnegative ray. It also proves P<0 for X>=7pi by complete root isolation and an exact endpoint sign.

On 11/2<=|x|^2<=7, complete Arb subdivision proves

f(x)<=M(2h(x)-h(x)^2).

There are 308 successful intervals, maximum depth 13. Rational Machin bounds enclose pi; centered polynomial Taylor bounds and monotone exponentials certify each interval. For |x|^2>=7 the physical function is negative while 0<h<=1, so the same coupled inequality follows automatically. The normalization H(0)>0 is checked exactly.

## PSD kernel, origin and density normalization

Outside the excluded ball, define

K(x,y)=(M/B)h(x)h(y)Q(xhat dot yhat).

Since Q has constant normalized harmonic coefficient b0=1 and all coefficients are nonnegative, the addition theorem gives PSD angular features. Set the origin degree-zero feature to -sqrt(MB), the nonorigin degree-zero feature to sqrt(M/B)h(x), and all positive-degree origin features to zero. Then

K(0,0)=MB, K(x,0)=-Mh(x), K(x,x)=Mh(x)^2.

The exact rank-one Gram relation (MB)(M/B)=M^2 is verified. A smooth radial cutoff chi, zero near 0 and one for |x|>=sqrt(11/2), gives a globally smooth rapidly decreasing extension: use sqrt(M/B)chi*h-sqrt(MB)(1-chi) in degree zero and chi*h times the positive-degree angular features. Required pair and triangle values are unchanged.

The complete cyclic sign established above and the coupled radial inequality verify the cyclic continuum hypotheses. This uses the cyclic version, not the literal individual-sign premise of Theorem 1.1. Dividing both f and K by H(0) gives

Delta_42 <= pi^21*(11/8)^21/21! * (P(0)+MB)/H(0).

The features are rapidly decreasing, so periodic sums converge absolutely; periodic approximation gives an unrestricted packing bound, not a lattice-only bound. For final rounding, the ball volume is pi^21*(11/8)^21/21!. A rational upper Machin bound prove the density is strictly below 0.000003677844262, independently of the Arb display enclosure.

## Discovery provenance

The angular LP was solved in the dimension-42 normalized harmonic basis, with sampled full cyclic constraints. Rationalization subtracts exactly 1/1000 from the constant coefficient before renormalizing it to 1. Only the subsequent complete exact/interval checks establish the new angular result.

The radial search used a scaled frozen44D degree-111 source polynomial plus a degree-31 linear combination of dimension-42 Laguerre Fourier eigenfunctions. Initial numerical radial outputs were invalid and had forbidden Fourier and physical-tail roots. Eight exact-root-guided repair rounds produced the frozen rational witness. The first angular candidate failed the cyclic proof; exact-critical-point separation produced the final accepted angular witness. M includes an explicit one-part-per-million upward guard. No optimizer status, sampled sign, or discarded tiny violation contributes to the proof.
