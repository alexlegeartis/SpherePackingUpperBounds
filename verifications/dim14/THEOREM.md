# Curved radial damping for the cyclic sphere-packing cone

This note concerns a new admissible kernel family. A numerical optimizer output is not a sphere-packing bound. The angular gates in `angular_verified.json` are rigorous; a density claim additionally requires the full radial and independent audits.

## Kernel and normalization

Let n=14 and minimum squared separation s=7/4. Put

h(u) = exp(-a(u/s-1)-b(u/s-1)^2), a=14, b=8, for u>=s.

Let Q(t)=sum_l c_l P_l^n(t), where the normalized Gegenbauer polynomials satisfy P_l^n(1)=1, c_l>=0, c_0=1, and B=Q(1). For M>=0 set

K(x,y)=(M/B) h(|x|^2)h(|y|^2)Q(<x,y>/(|x||y|)),
K(0,x)=K(x,0)=-M h(|x|^2), and K(0,0)=MB.

The degree-zero origin extension is the rank-one Gram form (M/B)(B,-h(x))(B,-h(y)); higher harmonic blocks are positive semidefinite by c_l>=0. For a smooth global extension, choose radial chi=0 near zero and chi=1 at radii at least sqrt(s), and eta=1 near zero and eta=0 at radii at least sqrt(s). Use degree-zero feature sqrt(M/B)*chi*h-sqrt(MB)*eta; multiply all higher angular features by chi. Quartic damping gives Schwartz features. All required origin and allowed-radius values are unchanged. The kernel need not be Fourier transformed.

Take f(x)=P(pi|x|^2)exp(-pi|x|^2), with nonnegative Fourier transform and H(0)=fhat(0)>0. The pair gate is

f(x) <= M(2h(|x|^2)-h(|x|^2)^2) for |x|^2>=s.

The cyclic gate is the sum of the three off-origin kernel terms of every triangle with all squared sides at least s. The standard periodic finite-sum argument then gives density at most

(pi*s/4)^(n/2)/Gamma(n/2+1) * (P(0)+MB)/H(0).

This is the cyclic finite-sum theorem, not a claim that each individual off-origin kernel term is nonpositive. Summing each rooted-star PSD form, pair terms are counted in their ordered-pair normalization, and every unordered triangle contributes twice its cyclic sum. The quartic radial damping is absolutely summable on separated point sets, so truncation/averaging introduces vanishing boundary errors. Only f uses Poisson summation.

## Reduction of all dilations to minimum side

Write sorted squared side ratios as x<=y<=z. Let q_0,q_1,q_2 be Q at the angles opposite z,y,x, respectively. After dividing by the positive factor h(sx)h(sy), the cyclic sum has sign

q_0 + r_1 q_1 + r_2 q_2,
where r_1=h(sz)/h(sy), r_2=h(sz)/h(sx).

Assume Q<=0 on [-1,tau] and Q>=0 on [tau,1], for some tau>=1/2. Then q_0<=0 and q_1>0 implies q_2>=0. Under a common squared-side dilation lambda>=1, both r_1 and r_2/r_1 decrease. Indeed, for v>=u>=1,

phi(lambda*v)-phi(lambda*u)
= lambda*(v-u)*(a-2b+b*lambda*(v+u)),

where phi(v)=a(v-1)+b(v-1)^2. Its lambda derivative is

(v-u)*(a-2b+2b*lambda*(v+u))>=0.

If q_1,q_2>=0, both positive weighted terms decrease. If q_1,q_2<=0, the total is automatically nonpositive. In the remaining case q_1<=0<=q_2, write the last two terms as r_1*(q_1+(r_2/r_1)*q_2). The bracket decreases; while positive its product with decreasing r_1 decreases, and after crossing zero it remains nonpositive. Thus the minimum-side=s cyclic inequality implies the inequality at every dilation. This also proves the lemma for any profile exp(-phi(u)) whose dilation differences phi(lambda*v)-phi(lambda*u) are nondecreasing.

## Minimum-side parameterization and certified domain

The side lengths are sqrt(s), u*sqrt(s), w*sqrt(s), with u>=1, w>=u. Write t for the angle between the first two sides. Then

w^2=1+u^2-2ut, -1<=t<=1/2, 2ut<=1.

The normalized cyclic expression is

Q(t) + r_1 Q((1-ut)/w) + r_2 Q((u-t)/w),

r_1=exp(-a*d-b*d*(2u^2-2+d)), d=1-2ut,
r_2=exp(-a*e-b*e^2), e=u^2-2ut.

`verify_angular.py` checks normalized Gegenbauer positivity and exact rational reconstruction. Exact Sturm sequences prove Q<=-1/20000 on [-1,1/2], exactly one real root in (1/2,1), and Q(1)>0. Thus the sign-transition assumption is a proved gate, not an inference from a graph.

The entire compact rectangle [-1,1/2] x [1,4] is covered by rigorous Arb enclosures. Rectangles are discarded only when their exact lower endpoints already force 2ut>1. Accepted rectangles use direct evaluation or a first-derivative mean-value enclosure. Rational bisection covers both endpoints; no numerical sign tolerance is used for acceptance. The source has independent corner counterexample checks before deep subdivision.

For u>=4, e>=u^2-1>=15, so r_2<=exp(-2010). If t>=-1/4, then w>=u gives (1-ut)/w<=1/u-t<=1/2, making the second Q term nonpositive. If t<=-1/4, then d>=3, r_1<=exp(-42). Since |Q(t)|<=B on [-1,1] by its nonnegative normalized Gegenbauer expansion, it suffices that

B*(exp(-42)+exp(-2010)) < 1/20000.

The verifier proves this with Arb. Hence no unbounded geometric domain is omitted.

## A useful obstruction to overly slow profiles

A noncompact inverse-power profile h(r^2)=r^(-2p) cannot satisfy the cyclic theorem with a nonzero Gegenbauer-positive Q and c_0>0. For collinear side lengths 1,r,r+1, divide by h(1)h(r^2). As r tends to infinity, the expression tends to Q(-1)+Q(1)=2*sum_even c_l>=2c_0>0, because h((r+1)^2)/h(r^2) tends to one and h((r+1)^2)/h(1) tends to zero. The same obstruction applies to any decaying radial profile with that consecutive-radius ratio tending to one. This rules out an apparently attractive algebraic-tail enlargement without an expensive numerical search.

## Remaining scope and reproduction

`screen.py` generated the angular candidate and initial radial pricing. Its numerical rows are discovery only. `verify_angular.py` rationalizes the stored source coefficients, applies the explicit 1/10000 constant guard, and checks the complete angular theorem. `search_radial.py`, `audit_radial.py`, and root-guided histories handle the radial side separately. All intermediate rejected candidates are retained. No family optimality or global sphere-packing optimality is asserted.
