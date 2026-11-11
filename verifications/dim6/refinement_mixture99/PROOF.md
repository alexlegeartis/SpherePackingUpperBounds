# Six-dimensional density certificate: proof and parameters

This frozen certificate proves unrestricted packing density **<0.409135399**. The finite-dimensional family and periodic argument are detailed below. All computational hypotheses pass full-domain verification; no inference rests only on numerical sampling.

## Frozen parameters

- Packing dimension n=6; minimum center separation 1; support radius R=3/2.
- Gaussian–Laguerre polynomial degree 31, radial Chebyshev degree 4, angular Gegenbauer degrees 0 through 12.
- Exact count radii 1.002005, 1.004023, 1.006053, 1.008096; corresponding average counts 79,80,81,82. Nonnegative exact weights are stored in witness.json.
- These counts follow by monotonicity from epsilon thresholds .002006,.004024,.006054,.008097 in [de Laat, de Oliveira Filho and Vallentin, Table 2](https://arxiv.org/html/1206.2608). Their half-open endpoint convention is preserved. The original spherical SDP certificates are external mathematical inputs.
- No Gaussian repair is needed: eta=zeta=0, with c=157/50 retained only for the general verifier. Exact rational Sturm checks prove the physical tail and Fourier sign globally. Both island lists are empty.
- The Fourier origin is exactly 500000000000001068179563/500000000000000000000000.
- The original positive-definite degree-zero block is saturated by the exact Schur operation. The remaining harmonic matrices and all off-origin function values stay unchanged.
- The density volume factor is pi^3/384. The exact saturated objective is in prerequisites.json and is independently reconstructed by verify.py.

The complete verification checks the pair interval, every compact cyclic triangle and every one-long-edge configuration, including support boundaries. The far/cyclic verifiers use exact Fraction geometry boxes, exact rational Bernstein initialization and outward-rounded float64 dyadic coefficient subdivision. The pair verifier uses 256-bit Arb Taylor bounds. FLINT/Arb, NumPy, correctly rounded binary64 arithmetic, and the stated analytic facts are trusted components.

The source-reviewed independent audit also reconstructs the full Fourier polynomial via a radial-Laplacian identity, checks all origin/guard arithmetic, and separately verifies the entire antipodal boundary with exact rational Bernstein coefficients. It is not a clean-room reimplementation of every continuum checker or a proof-assistant formalization.

The analytic framework is the existing cyclic three-point summation mechanism in [Cohn, de Laat and Salmon](https://arxiv.org/abs/2206.15373), augmented by the published average count inequalities. No novelty of this general theorem is asserted.

## Analytic family and density argument

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

## Discovery and exact source provenance

The frozen source is the exact rational convex combination of 1/100 of the
previous 0.410189987 source and 99/100 of a new joint-search checkpoint.
Discovery optimized the degree-31 radial polynomial and cyclic kernel jointly,
fixed an origin Schur relation while optimizing, lowered artificial PSD
regularization, and then froze the new radial polynomial for separation.
This improves discovery in the same analytic family; it is not a new
hierarchy theorem.

Before mixing, the new checkpoint's finite-decimal matrices were symmetrized
by exact averaging, and its isolated origin entry was replaced by the exact
Schur value plus 1/10^12. Every off-origin block already passed exact positive
definiteness, so no diagonal shift was needed. Those operations by themselves
do not certify the geometric inequalities. The new checkpoint failed a tiny
far-domain sign condition; only the final rational mixture is endorsed.

The final certificate is verified directly from its exact rational entries.
Its validity does not rely on trusting a solver, numerical separation, an
asserted convex-mixture result, or either endpoint as a complete certificate.
The serialized discovery history is provenance only. The original published
near-neighbor estimates remain external mathematical dependencies.
