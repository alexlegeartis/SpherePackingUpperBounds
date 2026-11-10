# Four-dimensional unrestricted packing-density certificate

The frozen witness proves density strictly below **0.635766243**, subject to the standard Fourier/positive-definite-kernel facts and the published neighbor-count bounds listed below. The computation has been independently audited within this research effort; it is not an external peer review or a proof-assistant formalization. No optimality or worldwide-priority claim is made.

## Exact witness and functions

Normalize minimum center separation to 1. Set R=333/200, and write X=pi*r^2. The 32 finite-decimal radial coefficients c_k in `witness.json` are interpreted as exact rationals. Define

f(x) = exp(-pi*|x|^2) sum_(k=0)^31 c_k L_k^1(2*pi*|x|^2)/(k+1).

For the Fourier convention with exp(-2*pi*i*x·y), the transform is the same expression with c_k replaced by (-1)^k c_k. This is the Gaussian–Laguerre Fourier eigenfunction identity in R^4. In particular H0=fhat(0)>0. The radial function is Schwartz, so Poisson summation applies.

The exact rational matrices Q_0,...,Q_16 are in `witness_psd_exact.json`. Their entries, not the original floating-point SDP matrices, define the kernel. For r in [1,R], put z=(2*r-R-1)/(R-1), and w(r)=(T_0(z),...,T_3(z)). For nonzero x in this annulus, set v_0(x)=(0,w(|x|)); at x=0 set v_0(0)=(1,0,0,0,0); otherwise set v_0(x)=0. For k>0 use w in the annulus and zero outside. Define

K(x,y) = v_0(x)^T Q_0 v_0(y) + sum_(k=1)^16 [w(|x|)^T Q_k w(|y|)] U_k(xhat·yhat)/(k+1),

with every k>0 term zero unless both arguments are nonzero and in the closed annulus. The degree-zero summand uses the same convention. The normalized Gegenbauer polynomial in four dimensions is U_k(t)/(k+1). Exact matrix positive definiteness and the spherical harmonic addition theorem imply that K is a positive-semidefinite kernel on all R^4. Continuity at the support boundary is not required for the finite kernel sums used below.

## External neighbor-count input

Let A_i be the average number of other packing centers at distances in [1,R_i). We use A_i<=B_i for

| R_i | B_i |
| --- | --- |
| 1.008096 | 24 |
| 1.017445 | 25 |
| 1.025977 | 26 |
| 1.036950 | 27 |

These follow by monotonicity from Table 2 of de Laat, de Oliveira Filho and Vallentin, [Upper bounds for packings of spheres of several radii](https://arxiv.org/pdf/1206.2608), whose corresponding epsilon thresholds are .008097, .017446, .025978, .036951. The epsilon-tangency graph is defined using distances [1,1+epsilon); its average-degree bounds are explicitly stated to be rigorous. The spherical SDP certificates underlying those published inputs are not reproduced in this package. Thus the replay is self-contained for the new witness, with these four cited results as external mathematical inputs.

The four lambda_i in `witness.json` are nonnegative exact finite-decimal rationals. The full local-pair inequality is

f(x)+2*K(x,0)+K(x,x) <= sum_i lambda_i 1_{|x|<R_i}, for |x|>=1.

The verifier proves it on [1,R] by Arb Taylor bounds, split at every R_i. At a threshold, the closed interval to its right verifies the smaller half-open allowance. Outside R, K and the allowances vanish and the exact tail sign supplies the inequality.

## Global cyclic condition

For every triangle with sides a,b,c>=1, the sum of the kernel at its three vertices is nonpositive:

C(a,b,c)=K(a,b,(a^2+b^2-c^2)/(2ab))+K(a,c,(a^2+c^2-b^2)/(2ac))+K(b,c,(b^2+c^2-a^2)/(2bc)) <=0.

Here the radius-angle notation is used and terms outside the kernel support are zero. If all sides are at most R, the verifier proves the polynomial numerator (abc)^16*C(a,b,c) negative over [1,R]^3. Every triple in this cube is a triangle because R<2. Symmetry permits restriction to a<=b<=c. The denominator is strictly positive.

If exactly one side exceeds R, only the opposite kernel term can remain. The far verifier proves K(r,s,t)<=0 on r,s in [1,R] with r^2+s^2-2rs*t>=R^2. If two or three sides exceed R, all terms vanish. Boundary cases are included in the closed-domain verifications.

## Density argument

Consider a periodic packing with N centers in a period cell of volume V. Enlarge the period lattice if necessary. Sum the PSD kernel over the finite displacement set consisting of zero and all supported neighbors around each center. After summing over centers modulo the period lattice, off-diagonal terms group into twice the cyclic sums of unordered triangles, hence contribute at most zero. Therefore

sum_(ordered distinct centers) [2*K(x-y,0)+K(x-y,x-y)] >= -N*K(0,0).

The pair inequality and the average count bounds now give

sum_(ordered distinct centers) f(x-y) <= N*K(0,0)+N*sum_i B_i*lambda_i.

Poisson summation and fhat>=0 give N^2*H0/V <= sum_(all ordered centers) f(x-y). Including the N diagonal terms yields

N/V <= [f(0)+K(0,0)+sum_i B_i*lambda_i]/H0.

Multiplication by the radius-1/2 ball volume pi^2/32 gives the claimed density bound. The supremum packing density is approximated by periodic packings: take arbitrarily large finite cube patches and periodically repeat them with a separation gap. The boundary-volume loss tends to zero. Thus the result is unrestricted, not lattice-only.

This is the cyclic three-point summation mechanism underlying [Cohn–de Laat–Salmon](https://arxiv.org/abs/2206.15373), augmented with the published local-count inequalities. No novelty of this general mechanism is asserted.

## Replay and rounding

`python verify.py` checks exact rational LDL, exact Sturm signs, the global pair inequality, both global kernel domains, and the final outward density. It fails if any condition is false. Bernstein coefficients are initialized exactly over the rationals, enclosed using Arb, and subdivided with outward IEEE-754 binary64 rounding via numpy.nextafter. Arb is used for domain exclusion. Standard correctly-rounded IEEE-754 arithmetic, NumPy and python-flint/Arb are trusted computational components; this is not a formally verified implementation.

The strict rational-polynomial checks have zero relevant roots. The pair proof uses 115 boxes; the far proof 34,009 boxes; the cyclic proof 83 boxes. The final Arb density enclosure is centered at 0.6357662426387588618572394093758579880966147479317909183832447245178820105673 with radius below 8e-77. The reported decimal .635766243 is rounded upward. The comparison uses that conservative decimal, not extra unreported digits.
