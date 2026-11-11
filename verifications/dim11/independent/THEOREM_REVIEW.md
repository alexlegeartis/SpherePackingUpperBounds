# Independent review of the curved-profile theorem

This audit uses the direct periodic cyclic argument, not the stronger individual-sign hypothesis in a printed three-point theorem. Let the minimum squared center separation be s>0 and define, for u>=s,

h(u)=exp[-a(u/s-1)-b(u/s-1)^2], with a,b>=0.

All four audited instances have a,b>0. Write Q=sum c_l P_l^n with normalized Gegenbauer polynomials P_l^n(1)=1, c_l>=0 and c_0=1, and put B=Q(1)>0. With M>0 define nonorigin K=(M/B)h(u)h(v)Q(t), cross K(0,x)=-Mh(u), and origin K(0,0)=MB. The constant harmonic is a rank-one Gram kernel with origin coordinate sqrt(MB) and nonorigin coordinate -sqrt(M/B)h(u). All other harmonic blocks are PSD by the addition theorem. Features below the packing separation may be set to zero; only the value at the isolated origin is special.

## Summation and density normalization

For a finite set of packing centers, sum the PSD quadratic form over each rooted star. Its expansion is

N*K(0,0) + sum_ordered_distinct_pairs [2K(x,0)+K(x,x)] + 2*sum_unordered_triangles C,

where C is the sum of the three kernel terms based at the triangle vertices. The factor two comes from interchanging the two nonorigin neighbors at each root. If C<=0 and f(x)+2K(x,0)+K(x,x)<=0 at every allowed distance, this expression is at most N*K(0,0)-sum_ordered_distinct_pairs f(x). It is nonnegative by PSD.

For a periodic packing with N centers per fundamental region of volume V, use the same expression per fundamental region. The sums converge absolutely: |Q|<=B by its nonnegative normalized Gegenbauer expansion, and h decays at least Gaussian (quartic in these cases), uniformly summably over separated sets. Thus the infinite rooted forms are limits of finite nonnegative forms, and triangle contributions can be regrouped exactly. Including diagonal f(0),

sum_{p in region} sum_{q in packing} f(q-p) <= N*(f(0)+MB).

Poisson summation for the polynomial-Gaussian f gives the opposite lower bound N^2*fhat(0)/V, because every Fourier coefficient fhat is nonnegative. Therefore center density N/V is at most (f(0)+MB)/fhat(0). Multiply by the ball volume at radius sqrt(s)/2. Arbitrary packing upper densities are bounded by periodic ones using large finite boxes and a separation-sized boundary gap before periodic repetition. Only f, not K, must be Schwartz for this argument.

Since Q(1)=B, K(x,x)=M*h(u)^2, so the precise pair condition is f(x)<=M*(2h-h^2). The physical polynomial is strictly negative from X=12, while h is in (0,1] at every allowed radius. Thus a complete compact pair proof through X=12 covers the whole infinite pair domain.

## Reduction of all triangles to minimum side

Sort squared side ratios x<=y<=z, all >=1. Let q_0,q_1,q_2 be Q at the angles opposite z,y,x. The largest angle has cosine <=1/2, hence q_0<=0. After division by h(sx)h(sy), the cyclic expression is

q_0 + r_1*q_1 + r_2*q_2,

where r_1=h(sz)/h(sy) and r_2=h(sz)/h(sx). The exact cap and one-transition polynomial gates establish Q<=0 below a single tau>1/2 and Q>=0 above it. Therefore q_1>0 implies q_2>=0; q_2<0 implies q_1<=0.

Under common squared-side dilation lambda>=1, for v>=u>=1 the exponent difference is

phi(lambda*v)-phi(lambda*u)=lambda*(v-u)*(a-2b+b*lambda*(v+u)).

Its derivative is (v-u)*(a-2b+2b*lambda*(v+u))>=0. Consequently r_1, r_2, and r_2/r_1 decrease. If q_1,q_2>=0 their contribution decreases. If both are nonpositive, the whole expression is nonpositive. In the remaining case q_1<=0<=q_2, write their contribution as r_1*(q_1+(r_2/r_1)*q_2). Its bracket decreases; while positive, the product decreases; once nonpositive it cannot cause a positive total. Thus the minimum-side=s inequality implies every dilation. This depends on the proved sign-transition gate, not on monotonicity of the whole polynomial.

## Complete minimum-side domain and unbounded tail

Write side lengths as sqrt(s), U*sqrt(s), W*sqrt(s), with U>=1, W>=U. Let T be the cosine between the first two sides. Then W^2=1+U^2-2UT, -1<=T<=1/2, and 2UT<=1. The three cosines are T, (1-UT)/W, (U-T)/W. The ratio exponents are

r_1=exp[-a*d-b*d*(2U^2-2+d)], d=1-2UT,
r_2=exp[-a*e-b*e^2], e=U^2-2UT.

The independent AD implementation covers the entire rectangle [-1,1/2] x [1,4], discarding only boxes with exact lower endpoints 2*T_min*U_min>1. It proves every retained point by interval bounds, not sampled evaluations. Although some retained boxes include unsorted triangles, every point still describes a real triangle. The identity W^2=1+U(U-1)+U(1-2T)>=1 justifies intersecting interval denominator enclosures with [1,infinity). Cosine interval clipping to [-1,1] is valid for the entire box. All eight partial derivatives have also been checked symbolically.

For U>=4 and the sorted domain, e>=U^2-1>=15, hence r_2<=exp[-15a-225b]. If T>=-1/4, the middle cosine is <=1/2: if 1-UT<0 it is already negative; otherwise W>=U gives (1-UT)/W<=1/U-T<=1/2. If T<=-1/4, then d>=3 and r_1<=exp[-3a]. The cap margin delta=1/20000 and |Q|<=B therefore imply the cyclic sum is negative whenever

B*(exp[-3a]+exp[-15a-225b])<delta.

The independent replay checks this strict inequality with rigorous Arb bounds for each instance. The exponents are dimension/profile specific: (24,3720), (30,3750), (30,3750), and (42,2010) for dimensions 11,12,13,14 respectively. Degenerate triangles and all boundary equalities are included.

## Exact volume

For n=2k+1, the ball volume is pi^k*s^(n/2)/(2^k*n!!). For n=2k it is pi^k*s^k/(4^k*k!). The independent verifier squares these positive quantities and uses rational alternating-series Machin bounds for pi, avoiding non-rigorous decimal rounding and correctly distinguishing odd and even dimensions.

Conclusion: no theorem, summation-factor, all-scale, tail, or normalization defect was found. The scripts certify explicit witnesses only, not optimality of any function family.
