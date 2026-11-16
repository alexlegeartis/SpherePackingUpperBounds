# Exact projective caps of arbitrary even degree

Let d>=3 and let X be N distinct unoriented lines, represented by unit vectors,
with |<x,y>|<=c<1 for distinct lines. Let P_k be the dimension-d normalized
Gegenbauer polynomial, so P_k(1)=1. Choose an integer r>=0 such that

    c^2 < (2r+1)/(d+2r).

Define mu_0=1 and mu_(2r)=(2r-1)!!/[d(d+2)...(d+2r-2)]. Then

    N <= floor((1-c^2)/(mu_(2r)*((2r+1)/(d+2r)-c^2))).

Proof. Use f(x)=x^(2r)(x^2-c^2). It is nonpositive for |x|<=c.
Write x^(2r)=sum_(j=0)^r a_(r,j) P_(2j)(x). Every coefficient is positive.
With lambda=(d-2)/2 the standard monomial expansion gives explicitly

    a_(r,j) = (2r)! Gamma(lambda) (2j+lambda) C_(2j)^lambda(1)
              / [2^(2r) (r-j)! Gamma(r+j+lambda+1)].

Taking the ratio of this expression at r+1 and r gives

    a_(r+1,j)/a_(r,j)
      = (r+1)(2r+1)/[(r+1-j)(d+2r+2j)].

Its denominator minus its value at j=0 is j(2-d)-2j^2<=0.
Consequently the ratio is smallest at j=0, where it equals
(2r+1)/(d+2r). The assumed strict inequality therefore makes every old
Gegenbauer coefficient of f positive; its new highest coefficient is positive
as well. Its constant coefficient is

    f_0 = mu_(2r)*((2r+1)/(d+2r)-c^2)>0.

Positive definiteness of each Gegenbauer kernel yields
f_0 N^2 <= sum_(x,y in X) f(<x,y>) <= N f(1).
Since f(1)=1-c^2, division proves the claimed bound. No numerical spherical-code
solver or approximate positivity test enters this argument.

For a centrally symmetric decomposition-point set, N counts antipodal pairs;
its number of points (and hence ordered source decompositions) is at most 2N.
Take the minimum over all admissible r if desired. Cases r=0,1,2 are the relative,
quartic and sextic caps. Their pair bounds respectively are

    d(1-c^2)/(1-d c^2),
    d(d+2)(1-c^2)/(3-(d+2)c^2),
    d(d+2)(d+4)(1-c^2)/(3*(5-(d+4)c^2)).

For every c<1 there exists an admissible r: it suffices that
r>(d*c^2-1)/(2*(1-c^2)). This ensures finite exact multiplicity caps in every
fixed dimension, though it does not ensure the resulting density cut is strong.
