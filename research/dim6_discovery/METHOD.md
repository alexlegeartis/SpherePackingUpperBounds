# Dimension-generic Chebyshev/count cyclic family

The following analytic family is valid. Individual numerical candidates are not valid until every required condition is certified.

Normalize minimum center separation to one. Let 1<R<2, let w(r) contain Chebyshev polynomials T_0,...,T_p in z=(2r-R-1)/(R-1), and set all radial features to zero outside [1,R]. Add an isolated origin coordinate in harmonic degree zero. Let Q_0 have size p+2 and Q_l (l>0) size p+1. Exact positive semidefiniteness of every Q_l, together with the spherical harmonic addition theorem, gives a globally positive-semidefinite kernel K. Its angular functions are the normalized Gegenbauer polynomials in the packing dimension n.

For the Fourier convention exp(-2 pi i x dot xi), set X=pi |x|^2 and

f(x)=P(X) exp(-X),
P(X)=sum_k c_k L_k^(n/2-1)(2X)/L_k^(n/2-1)(0).

The Fourier polynomial H is obtained by multiplying c_k by (-1)^k.

For nonnegative rational weights lambda_i and published valid average counts B_i within distances [1,r_i), require the complete pair inequality

f(x)+2K(x,0)+K(x,x) <= sum_i lambda_i 1_(|x|<r_i),  |x|>=1.

At every threshold the smaller right-hand allowance must hold. For every triangle with side lengths at least one, require that the sum C of K at its three vertices be nonpositive. The compact triangle region is [1,R]^3, all of whose points satisfy the triangle inequality when R<2. With exactly one long edge, the only surviving kernel must be nonpositive; with two long edges all kernels vanish. The pair condition outside support is the global physical tail sign. Require H>=0 globally and H(0)>0.

Periodic finite PSD summation, grouping each unordered triangle with multiplicity two, then Poisson summation yields packing density at most

vol(B_n(1/2)) * [P(0)+Q_0(0,0)+sum_i lambda_i B_i]/H(0).

The passage from periodic to unrestricted packings is the standard large-patch periodic approximation. This is the existing cyclic three-point mechanism with published count cuts, not a new general hierarchy theorem.

## Exact origin saturation

Write Q_0=[[a,b^T],[b,B]], with B positive definite. If alpha=a-b^T B^-1 b >=0, replace a by a-alpha. The resulting block is exactly PSD, through the congruence [[b^T B^-1],[I]] B [[B^-1 b,I]]. Every nonzero pair, compact triangle, and one-long-edge expression is unchanged. Only K(0,0) decreases by alpha. Its density cost reduction is vol(B_n(1/2))*alpha/H(0).

## Reciprocal Gaussian repair

If the source polynomial has forbidden distant sign islands, add

g(x)=eta(c-X) exp(-X/4)+zeta(c-X) exp(-4X),

with eta,zeta>=0 and n/8<c<pi. Its Fourier transform is

eta 2^n(c-2n+16Y) exp(-4Y)
 + zeta 2^(-n)(c-n/8+Y/16) exp(-Y/4),

where Y=pi|xi|^2. Both physical corrections are nonpositive for |x|>=1, so a previously verified pair inequality remains valid. The second Fourier correction is positive everywhere; the first is negative only for Y<(2n-c)/16.

Exact rational Sturm complements and outward interval bounds must prove both global signs. The changed Fourier normalization is

H0'=H(0)-eta 2^n(2n-c)+zeta 2^(-n)(c-n/8).

The changed numerator is P(0)+K(0,0)+sum lambda_i B_i+c(eta+zeta). Origin saturation must divide alpha by H0', not by the old H(0). Guarding a radial function does not repair a positive pair, cyclic, or far kernel value; those conditions need their own complete verification.

## Polynomial certification

The compact cyclic expression is rational in side lengths, with positive denominator (abc)^K for maximum angular degree K. Multiplying the complete cyclic sum by this denominator yields a rational polynomial. Exact Bernstein coefficient initialization and outward dyadic subdivision establish the complete sign. The one-long-edge check works in radii and angle and excludes only boxes proved geometrically impossible. The pair verifier uses exact rational source coefficients and 256-bit Arb Taylor enclosures, splitting at every count threshold.
