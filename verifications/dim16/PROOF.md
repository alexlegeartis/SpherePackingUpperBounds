# Combine a widened spherical count with a radial cyclic kernel

This document states the linkage between the separately checked components.
Final package freezing and independent radial audit are tracked separately.

## 1. The positive-width count

The exact local spherical certificate has cap c=5003/10000 in dimension16
and real objective less than7803. Its Schoenberg coefficients and all
reweighted local Gram blocks are nonnegative/PSD. Its global pair and
three-variable cyclic signs have been checked on their complete domains.
The usual finite-sum argument gives N<=7802.215167128724, hence the
integer code cardinality satisfies N<=7802.

Fix minimum squared separation s=2 and squared shell threshold
tau=10000/4997. Around a packing center, take all other centers with
s<=r^2<tau, and radially normalize them to unit length. If their squared
radii are u,v and their squared mutual separation is w>=s, the inner
product is at most (u+v-s)/(2 sqrt(u v)). This expression increases in
both u and v on [s,tau], since tau<2s and its derivative in u has the
sign of u-v+s. Thus its maximum is 1-s/(2tau)=5003/10000.
There are consequently at most7802 neighbors in this half-open shell.

This is a pointwise count here, not merely an average-degree input. At the
outer threshold itself its own pair allowance is absent. The inclusive
spherical cap proof makes this convention conservative.

## 2. The radial cyclic kernel

Let Q have nonnegative dimension16 Schoenberg coefficients and constant
coefficient one, with B=Q(1). For r^2>=s define
h(r)=exp(-phi(r^2/s)), where phi is the frozen nonnegative linear-hinge
and logarithmic-hinge combination. The origin is a separate feature:

    K(0,0)=C,
    K(0,x)=K(x,0)=-(C/B) h(|x|),
    K(x,y)=(C/B^2) h(|x|)h(|y|)Q(cos(x,y)).

Here C>=0. The degree-zero origin block is
C*(1,-1/B)^T*(1,-1/B), so is PSD. Every higher harmonic has a
nonnegative coefficient. Thus each rooted local quadratic form is
nonnegative. Values for unused radii below sqrt(s) may be extended with
h=1; continuity at the isolated origin is not needed in the finite-sum
argument, and K is never Fourier transformed.

The exact geometry proof verifies the sum of the three rooted off-diagonal
kernel terms at every packing triangle is nonpositive. It uses the full
signed cyclic condition, not three separate pointwise sign inequalities.
All scales and the infinite tail are included. Since h has a Gaussian
outer tail, all lattice-periodic local sums converge absolutely.

## 3. The radial function and count allowance

The radial function is f(x)=exp(-pi|x|^2)P(pi|x|^2), with its exact
dimension16 Fourier partner exp(-pi|y|^2)H(pi|y|^2). The certificate
proves H>=0 on the whole nonnegative ray and H(0)>0.
The complete pair condition is

    f(x) - (C/B)(2h(|x|)-h(|x|)^2)
         <= lambda * 1_{s<=|x|^2<tau},

where lambda>=0. The verifier handles the smaller allowance at equality
|x|^2=tau. It also verifies every profile breakpoint and the entire
remaining infinite tail. Physical negativity by itself is not assumed
at the contact radius; only the displayed pair condition is required.

## 4. Density bound

For a periodic packing, Poisson summation applied only to the Schwartz
function f gives the standard lower contribution rho*H(0) per center.
Add the nonnegative rooted kernel forms. The pair terms are at most
lambda times the ordered shell degree, whose average is at most7802.
The triangle terms are nonpositive by the cyclic certificate. Therefore

    rho <= [P(0)+C+7802*lambda]/H(0).

Multiplying by the volume of a ball of radius sqrt(s)/2 gives

    density <= (pi^8/8!)*(s/4)^8
               * [P(0)+C+7802*lambda]/H(0).

The same upper bound applies to arbitrary packings by the standard
periodic approximation argument. The volume factor at s=2 is
pi^8/10321920, not the unit-ball volume or the Barnes-Wall density.

## 5. Repair and normalization

The numerical radial output was normalized exactly, then repaired once
by adding epsilon=1/1000000 to its normalized degree-one Laguerre
coefficient. In dimension16 this adds
epsilon*(1-X/4) to P(X) and epsilon*(X/4-1) to H(X).
Thus H(0)=1-epsilon exactly before rescaling. In the frozen final export,
P, H, C and lambda are all divided by the same exact1-epsilon, giving
H(0)=1. This is equivalent to dividing the unscaled objective by1-epsilon;
C and lambda are not silently left at an inconsistent normalization. All
radial signs and the full pair inequality were reverified after the repair.

The strict outward-rounded primary result with U=7802 is

    density16 < 0.024920988.

The underlying Arb enclosure begins0.0249209879146585408. This is a small
improvement over the pinned published0.02492121, not a large closure of
the gap to the achieved Barnes-Wall density pi^8/645120.

No unsolved optimizer status, finite sampling, or source paper's numerical
table is used as proof of any new inequality. Independent verification
status and exact file manifests must accompany any eventual release.
