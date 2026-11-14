# Proof of the 46D upper bound

The conclusion is

    Delta_46 < 0.0000008756454737.

The underlying theorem is the existing cyclic Cohn–de Laat–Salmon continuum
framework ([arXiv:2206.15373](https://arxiv.org/abs/2206.15373), Section 2.2,
equation (2.4), and Theorem 2.4), together with its periodic PSD summation
argument. The cyclic sum condition is essential; the individual-sign premise
of the printed Theorem 1.1 is not substituted for it. This certificate has no
external neighbor-count input.

## Exact functions and normalization

Let s=6, kappa=17/8, X=pi*|x|^2, and h(x)=exp[-kappa*(X-6*pi)].
`certificate.json` gives rational degree-111 polynomials P and H and a positive
rational M. Set f(x)=P(X)*exp(-X). Both the direct coefficient-sum formula
and an independent Laguerre three-term recurrence prove exactly that H is
the dimension-46 Fourier polynomial of f. In particular, fhat(0)=H(0)>0.

The angular polynomial Q has degree 19, normalized dimension-46 Gegenbauer
coefficients b_l>=0, b_0=1 and B=Q(1)=sum(b_l)>0. Primary closed-form and
independent recurrence/moment calculations reconstruct these identities.
For allowed nonzero vectors put

    K(x,y) = (M/B)*h(x)*h(y)*Q(xhat dot yhat).

The addition theorem gives positive semidefinite nonorigin angular features.
In degree zero choose the origin feature -sqrt(M*B) and the nonorigin feature
sqrt(M/B)*h(x). Positive-degree origin features vanish. Thus

    K(0,0)=M*B,  K(0,x)=-M*h(x),  K(x,x)=M*h(x)^2.

The exact identity (M*B)*(M/B)=M^2 proves the rank-one Gram condition. A smooth
cutoff chi, zero near the origin and one for |x|>=sqrt(6), extends the features
globally: use sqrt(M/B)*chi*h-sqrt(M*B)*(1-chi) in degree zero and chi*h times
the positive-degree angular features. This preserves every required value
and makes the origin extension smooth. Features decay rapidly at infinity.

## Complete cyclic condition

Both angular implementations prove Q<-1/100000 on [-1,1/2], Q<0 on
[-1,5001/10000], and precisely one subsequent root before 1, using exact
rational Sturm chains. Nonnegative harmonic coefficients imply |Q|<=B.

For sorted normalized triangle coordinates 0<=t<=1/2, 1<=u<=5, 2*t*u<=1,
define w=sqrt(1+u^2-2*t*u), z1=(1-u*t)/w, z2=(u-t)/w. Complete interval
subdivision proves

    Q(t)+exp[-40*(1-2*t*u)]*Q(z1)
        +exp[-40*(u^2-2*t*u)]*Q(z2) < 0.

The primary check uses the inherited natural interval radicand. The independent
check uses exact monotone rational radicand bounds and nonnegative sorted-gap
bounds. It certifies 4212 rectangles and discards only
rectangles with 2*t*u>1. Subdivision is a full-domain proof.

The single-sign-transition scale lemma extends the result to every larger
shortest side and every alpha>=40. The actual parameter is
alpha=6*pi*kappa=51*pi/4>40. To see the lemma, write the secondary terms as
A*exp(-p*T)+C*exp(-q*T), 0<=p<=q and T>=1. The ordered angles give either two
nonnegative terms, two nonpositive terms, or A<=0<=C. In the mixed case,
factoring exp[-p*(T-1)] shows the sum is at most the maximum of its value at
T=1 and zero. Since Q(t)<0, each case preserves strict negativity.

For t<=0 the secondary exponential gaps are at least 1. For u>=5 a positive
first secondary contribution forces its gap above 4.001; the second gap is
at least 24. Exact arithmetic checks 40^30/30!>2*B*100000, bounding these
positive tails below the strict negative cap margin. These analytic arguments
cover all scales and the infinite triangle tail; the compact rectangle proof
alone is not the global theorem.

## Global radial and pair gates

Complete certified complex-root isolation accounts for all 111 roots of each
polynomial and excludes real roots of H on [0,infinity) and of P on
[40,infinity). Exact H(0)>0 and P(40)<0 establish the ray signs. Both
implementations use FLINT/Arb root isolation; they do not constitute separate
root-isolation backends.

The required pair inequality is

    P(X)*exp(-X) <= M*(2*h(X)-h(X)^2),  X>=6*pi.

It is equivalent, after division by positive h, to

    P(X)*exp[(9/8)*X-(51/4)*pi]
        <= M*(2-exp[-(17/8)*(X-6*pi)]).

The primary proof uses centered polynomial Taylor upper bounds and monotone
exponentials on an enclosing rational interval [6*pi_lower,40]. It accepts
9927 intervals, maximum depth 25.
The independent proof encloses the complete signed expression and explicitly
checks that its 9967 accepted leaves tile that
interval, maximum depth 25. Machin bounds obtained
from alternating arctangent sums enclose pi. Beyond 40, P<0 and 0<h<=1 make
the pair inequality automatic.

`tail_escape.json` gives an exactly positive P value at a rational point
beyond 7*pi. Thus the earlier stronger physical-negativity-from-7*pi template
cannot certify this witness; the extended pair proof is necessary. The
positive value is compatible with the actual coupled inequality.

## Density conclusion and strict rounding

The cyclic periodic summation and PSD hypotheses bound center density by
(P(0)+M*B)/H(0). Absolute convergence follows from rapid decay. Periodic
approximation extends the bound to unrestricted packings. Multiplying by
the radius-sqrt(6)/2 ball volume gives

    Delta_46 <= (3*pi/2)^(46/2)/Gamma(46/2+1)
                 * (P(0)+M*B)/H(0).

For even n the upper Machin bound for pi directly gives rational strict
rounding. For odd n=2m+1 the squared volume is
4*pi^(n-1)*6^n*((m+1)!)^2/((n+1)!)^2, so rational comparison of squared
positive quantities proves the same conclusion without relying on displayed
floating digits. The independent density enclosure is

    [8.756454736982146898517898346872551143894758022735253364759971983646334842359551627953120322406796137838004170727886344471917134092380775344247377150070852765174934082855413056597306941486389407814027036027970674765601287884367063e-7 +/- 7.13e-236]

Both strict-rounding checks prove the stated decimal upper bound. Numerical
LP objectives, dense samples and the search checkpoint's solver status are
not proof inputs. This is a computer-assisted certificate, without a claim
of optimality, formal proof-assistant verification or external peer review.
