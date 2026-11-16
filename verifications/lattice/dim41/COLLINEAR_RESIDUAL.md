# Exact collinear residual bridge

Under the hypotheses of CYCLIC_LATTICE_BRIDGE.md, partition the ordered distinct-nonzero lattice pairs into the set

 L={(v,-v),(v,2v),(2v,v):v in Lambda, v!=0}

and its complement. These three parametrizations are disjoint as ordered-pair sets. L is invariant under the cyclic transformations (x,y)->(-x,y-x)->(-y,x-y), hence its complement is also invariant. The complement's total K-sum is nonpositive by the certified cyclic inequality.

By symmetry of K, the sum over L is exactly

 sum_{v!=0} R(v), R(v)=K(v,-v)+2K(v,2v).

There is no extra factor two. The second copy represents the distinct ordered pair (2v,v). Both signs of v are already present in the lattice sum.

PSD and absolute convergence give

 C+sum_{v!=0}[D(v)+R(v)] >=0.

Consequently any Fourier-positive radial f and nonnegative compact envelope G satisfying f+D+R<=G on all |v|>=1 yield

 rho*fhat(0)<=f(0)+C+2(2^n-1)kappa(G).

The profile/antipodal counting factor remains exactly the same as before.

For the certified radial-harmonic kernel,

 R(t)=C/B^2 [Q(-1) h(t)^2+2B h(t)h(4t)],
 D(t)+R(t)=C/B [(1+Q(-1)/B)h(t)^2-2h(t)+2h(t)h(4t)].

For the recovered38D kernel, a=1+Q(-1)/B is between0 and1, and h decreases from1 on t>=1. Therefore a h^2-2h increases with t, while 2h(t)h(4t) decreases. On a subinterval [u,v], a rigorous upper bound is obtained by evaluating the main term at v and the positive cross term at u, with outward arithmetic.

The ratio h(4t)/h(t)=4^96 exp(-360t) is decreasing. Proving Q(-1)/B+2*4^96 exp(-360)<0 proves R(t)<0 for the entire ray t>=1. Thus D+R<=D<=0, so a separately certified negative polynomial tail for f is sufficient beyond the radial support.

Implementation audit: search.py uses mult=2(2^n-1), radial column mult/B times the displayed bracket, physical C as a nonnegative variable with objective coefficient1, and scaled profile budgets. These factors agree with this theorem. Current sampled Fourier negativity invalidates the pilot as a certified bound until repaired and replayed.
