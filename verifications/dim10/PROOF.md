# Reconstructed 10D log-hinge cyclic bound

This package is a fresh reconstruction, not recovery of byte-identical lost files.
It proves density < 0.143330282 for arbitrary sphere packings in dimension 10.
Its pinned comparison is the displayed 0.14340271 in arXiv:2512.04829v2,
Table 1, not an independently replayed comparator or worldwide-priority claim.

Set the minimum squared center separation to s=5/4. The exact angular
polynomial Q has nonnegative normalized dimension-10 Gegenbauer coefficients,
constant Gegenbauer coefficient 1, and B=Q(1). Its coefficients are frozen in
geometry/angular.json. The exact Sturm checks prove Q<-1/20000 on [-1,1/2]
and a single negative-to-positive sign transition on (1/2,1).

For v>=1 let phi(v) be the nonnegative linear and logarithmic hinge combination
in geometry/fit.json, using linear knots 1.02,1.05,1.1,1.2,1.4,2 and logarithmic
knots 1.2,1.4,1.5,1.6,1.65,1.75. Define h(r)=exp(-phi(r²/s)). All decimal
parameters are interpreted as exact rationals. The increasing function
v phi'(v), including nonnegative jumps at knots, makes
h(sqrt(lambda*z*s))/h(sqrt(lambda*y*s)) nonincreasing in lambda for z>=y>=1.

For sorted triangle side lengths 1<=u<=w, w²=1+u²-2ut, the normalized
cyclic expression is

    Q(t) + exp(phi(u²)-phi(w²))*Q((1-ut)/w)
         + exp(-phi(w²))*Q((u-t)/w).

The interval verifier covers [-1,1/2] x [1,4], discarding only boxes wholly
outside 2tu<=1. It proves all remaining boxes strictly negative, with exact
partition area 9/2. A centered mean-value interval bound uses one-sided
hinge derivatives enclosed together; local Lipschitz continuity justifies it.
For u>=4 the exact final linear slope gamma>=7 proves the infinite tail:
when t>=-1/4 only the last angular value may be positive and its contribution
is at most B exp(-14 gamma); when t<=-1/4 both remaining positive parts total
at most B(exp(-3 gamma)+exp(-14 gamma)). Both are <1/20000.
The ordered angular values and the one-sign-transition property propagate
the signed cyclic inequality to all larger common scales: if the two later
Q values are nonnegative their multipliers decrease; if the middle Q is
negative factor its multiplier and use the decreasing ratio of the last
multiplier to it; if both are negative no further argument is needed.

Use the rooted PSD kernel

    K(0,0)=C,
    K(0,x)=-(C/B)*h(|x|),
    K(x,y)=(C/B²)*h(|x|)*h(|y|)*Q(cos(x,y)).

The degree-zero block is C*(1,-1/B)^T*(1,-1/B). All higher harmonics are PSD
by the addition theorem. The complete triangle contribution is nonpositive
by the cyclic proof. This is a cyclic finite-sum argument, not a claim that
every single off-diagonal angular term is nonpositive.

The Schwartz radial function f(x)=exp(-pi|x|²)P(pi|x|²) has Fourier partner
exp(-pi|y|²)H(pi|y|²). The radial verifier independently decomposes P into
normalized Laguerre polynomials L_k^4(2X)/L_k^4(0) and verifies the exact
(-1)^k Fourier transform. H(0)=1 and H>0 on [0,infinity), proved with
complete exact-rational complex-root isolation and root-count checks.
P<0 for X>=7. Exact interval Taylor bounds prove

    exp(-X) P(X) <= (C/B)*(2h-h²)

for v=X/(pi*s) in [1,2], splitting at every profile knot. Since 2*pi*s>7,
physical negativity covers the full remaining ray. No external count
inequality or pointwise kissing bound is used.

For a periodic packing of center density rho, Poisson summation applied to f
and the nonnegative rooted forms yield rho <= P(0)+C: the origin contributes
P(0)+C, pair terms are nonpositive by the pair condition, and triangle terms
are nonpositive by the cyclic condition. Gaussian tail decay makes all sums
absolutely convergent. Standard periodic approximation extends the bound to
arbitrary packings. Multiplying by ball volume gives

    density <= (pi^5/120)*(5/16)^5*(P(0)+C)
            < 0.143330282.

The computed enclosure begins 0.14333028156576795078. The coefficients were
normalized, repaired with epsilon=1/1000000 in the normalized degree-one
Laguerre direction, renormalized, and C enlarged by 1000001/1000000. All
claimed inequalities were rechecked on the final coefficients; a numerical
optimizer status is not proof. The analytic ingredients are the standard
Gegenbauer addition theorem, Poisson summation and periodic approximation;
they are not formalized in this repository. The continuum framework follows
Cohn, de Laat and Salmon, https://arxiv.org/abs/2206.15373.
