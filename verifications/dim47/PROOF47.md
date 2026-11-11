# Dimension 47 by restriction and angular-origin saturation

## Result

The inherited 48-dimensional proof and the exact transfer checks prove the strict unrestricted sphere-packing density upper bound

**Delta_47 < 0.000000616879891.**

The exact bound has an Arb enclosure centered at

6.1687989093412065475544357585842309681792844038065763031987165432183814511319767474205840514692697459229260992961213371795782539452384696219647263669872e-7,

with enclosure radius below 6.5e-159 at 512-bit precision. This improves the pinned published comparator 0.0000006490757338 by approximately 4.96026%. The comparator is Cohn's indexed table citing Afkhami-Jeddi et al., arXiv:2006.02560v2. Direct access to the table returned 502 on 2026-11-11, while the search index returned a four-month-old crawl. This comparison does not establish worldwide priority.

This is a computer-assisted upper-bound proof, not a packing construction, proof of optimality, external peer review, or proof-assistant formalization. The mathematical density theorem is the existing Cohn–de Laat–Salmon cyclic continuum framework; see the inherited source proof and arXiv:2206.15373.

## Frozen source

All input coefficients are inherited, without numerical reoptimization, from the 48D certificate in alexlegeartis/SpherePackingUpperBounds at commit e756587e8fdf7827f085d6d0372cac267b110bd1. The directory `continuum_three_point` contains a byte-for-byte copy of the source 48D package. The two defining inputs are:

- `adaptive_cyclic40_radial111.json`: rational polynomials P and H of degree 111, with f_48(x)=P(pi|x|^2) exp(-pi|x|^2), Fourier transform H(pi|xi|^2) exp(-pi|xi|^2).
- `cyclic_angular40_exact.json`: rational degree-19 polynomial Q, B=Q(1)>0, initially expanded in normalized dimension-48 Gegenbauer polynomials.

Use the **unnormalized** coefficients from these inputs. Define the positive rational

M = neighbor_cost / B * (1000001/1000000).

The inherited source verifies H>=0 everywhere on the nonnegative ray; the exact Fourier identity; and, for all |x|^2>=6,

f_48(x) <= M (2h(x)-h(x)^2),

where h(x)=exp(-(17/8)pi(|x|^2-6)). Its proof includes the complete compact coupling interval, the infinite physical tail, the complete cyclic triangle domain, and analytic tails. Its angular/cyclic argument is independent of the ambient dimension once Q is fixed.

Fresh replay of the three inherited independent verifiers passed in this package: angular Sturm checks; radial identity, global signs and 33 coupling intervals; and 11,932 compact cyclic rectangles with 445 infeasible rectangles excluded, plus analytic tails. The cyclic proof uses alpha=51pi/4>40.

## Restriction preserves Fourier positivity

Restrict the physical function to R^47: f_47(x)=f_48(x,0). Since the source is Schwartz, Fourier inversion and Fubini give, with the 2pi Fourier convention,

hat f_47(xi) = integral over R of hat f_48(xi,eta) d eta >= 0.

Its origin value is the positive rational

N_47 = sum from k=0 to 111 of p_k (47/2)_k,

where (a)_k denotes the rising factorial and p_k is the coefficient of X^k in P(X). The exact computation checks N_47>0. It also reconstructs all 112 Fourier coefficients two ways: directly using the dimension-47 Gaussian Fourier identity, and by integrating H(pi|xi|^2+pi eta^2) exp(-pi eta^2). They agree coefficient by coefficient. The latter integral identity carries the source's full global Fourier positivity to dimension 47 without replacing it by a sampled sign test.

The physical function itself is unchanged as a function of radius. Thus its source coupling inequality still holds at the same minimum squared separation 6.

## The essential origin improvement

Simply restricting the entire old kernel would give a bound near 7.82168e-7 and would not beat the comparator. The following exact positive-semidefinite reconstruction is essential.

Expand the SAME angular polynomial in normalized dimension-47 Gegenbauer polynomials:

Q(t) = sum_l b_l P_l^(47)(t).

Exact triangular coefficient conversion gives every b_l>=0 and

b_0 = 1.3710462822361538... > 1.

The verifier independently computes this constant from spherical moments:

b_0 = sum_j [t^(2j)]Q(t) (1/2)_j / (47/2)_j.

For both nonzero arguments outside the excluded ball, keep precisely the old kernel:

K(x,y) = (M/B) h(x)h(y) Q(xhat dot yhat).

Its degree-zero feature is sqrt(M b_0/B) h(x). Choose its value at the origin to be -sqrt(M B/b_0), and let all positive-degree angular features vanish at the origin. Therefore

K(x,0) = -M h(x),
K(x,x) = M h(x)^2,
K(0,0) = M B / b_0.

All required nonorigin and cross values are unchanged, while the origin cost decreases from MB to MB/b_0. The rank-one scalar Gram identity

(MB/b_0) (M b_0/B) - M^2 = 0

is checked over the rationals. Nonnegative Gegenbauer coefficients and the addition theorem prove positive semidefiniteness. This is an exact origin Schur saturation, not a rounded singular matrix or arbitrary subtraction.

To obtain a globally smooth rapidly decreasing kernel, choose a smooth radial cutoff chi that is zero near 0 and equal to 1 for |x|>=sqrt(6). The degree-zero feature is

sqrt(M b_0/B) chi(x)h(x) - sqrt(M B/b_0)(1-chi(x)).

For l>=1 multiply the spherical harmonic features by chi(x)h(x). All are smooth near the origin and rapidly decreasing. No value in a required pair or triangle configuration is changed.

## Density theorem and normalization

The source's cyclic triangle inequality holds for exactly the same nonorigin kernel values, so it remains valid. Its coupled pair condition uses the same K(x,0), K(0,x), and K(x,x), so it too is unchanged. The only modified theorem quantity is K(0,0). Divide both f_47 and K by N_47. Then hat f_47(0)=1, all hypotheses are satisfied, and the cyclic continuum density theorem gives

Delta_47 <= pi^(47/2) (3/2)^(47/2) / Gamma(49/2) * (P(0)+MB/b_0)/N_47.

The verifier evaluates this expression with Arb intervals and proves it strictly smaller than the displayed decimal bound. The exact coefficients and expressions, not floating approximations, define the witness.

## General transfer and limitations

For every integer n<=48 for which the transferred coefficients remain positive semidefinite, replace 47 by n in the formulas. In fact the exact conversion here checks b_l>=0 for every n=25,...,48. Only dimension 47 beats its pinned comparator in the direct transfer screen; dimension 46 is about 2.02% weaker, and dimensions 25–45 are weaker still. This does not prove that the broader method cannot improve those dimensions. It only describes the specific inherited witness.

Two additional numerical one-parameter screens were tested: multiplication by an extra physical Gaussian, and Gaussian-weighted integration over discarded coordinates. They did not improve the direct restriction outcome. These optimizer screens are explicitly not infeasibility proofs and are not used in the certified 47D result.
