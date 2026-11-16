# From the trace theorem to the density certificate

All vector radii below are squared norms. Normalize the lattice minimum to1. Put n=18, tau=20 and h(x)=p(tau|x|²)exp(−tau|x|²/2). The polynomial is expanded in generalized Laguerre polynomials of parameter8, normalized at zero. Fourier transformation changes coefficient k by(−1)^k; the resulting polynomial q is certified strictly positive on[0,infinity). Consequently h-hat is nonnegative and h-hat(0)=(2pi/tau)^9 q(0).

The verifier proves p<0 on[45,infinity), so h<0 at squared norms≥9/4. It verifies a rational constant majorant g_b on every partition interval[E_b,E_(b+1)]. To avoid transcendental arithmetic, write a=tau E_b and u=tau(t−E_b): the polynomial inequality p(a+u)≤Gamma_b(1+u/2+u²/8) implies h≤Gamma_b exp(−a/2). Replacing exp(a/2) by its positive128-term Taylor lower bound weakens this to a rational bound. Long-interval heights are raised if necessary to be nonincreasing.

## Class cardinality constraints

For distinct antipodal pairs in the same nonzero mod2 class, both half-sums are nonzero lattice vectors. Thus |v·w|≤(|v|²+|w|²−4)/2. Pairs of squared norm≤s have normalized absolute inner products≤c=1−2/s. The relative bound on antipodal lines gives

    B(s)=floor[n(1−c²)/(1−n c²)].

There is at most one pair of norm below2 in a class. If one such pair has norm≤a<2, no other pair of norm≤s exists when a+s≤4. Otherwise its projected companions obey the corresponding(n−1)-dimensional relative bound with

    b=(a+s−4)²/(4as), c'=(c+b)/(1−b).

The allowed number of companions is the minimum of B(s)−1 and that projected relative bound when its denominator is positive. Upper endpoints a,s give valid weak constraints throughout each bin. Backward minima of these cumulative caps remain valid.

## Complete first18 profile disjunction

Assume rho>32768, so TRACE_THEOREM.md supplies trace excess≥1/25 for every18 same-class pairs. In an empty-short class put t=18, shortdev=shortheight=0. In a class with a short pair in bin j put t=17, shortdev=E_(j+1)−2, shortheight=g_j. Let x_b be the nonnegative number of long pairs in bin b, and a_b=E_(b+1)−2.

The cumulative bounds are sum_(b≤k)x_b≤cap_k. The light branch has total long count≤t−1. For the heavy branch let h be the bin of the t-th long pair. Then the following relaxed constraints are valid:

    sum_(b<h)x_b≤t,
    sum_(b≤h)x_b≥t,
    sum_(b<h)(a_h−a_b)x_b≤t a_h+shortdev−1/25.

The last inequality bounds the excess of the first18 actual pairs above by the bin endpoints, then uses the trace lower bound. The before-h count is deliberately relaxed from t−1 to t. Cutoffs with cap_h<t or negative final right-hand side are impossible. Every other cutoff is represented by a rational dual certificate.

For nonincreasing g, cumulative-cap greedy filling gives a valid bound on the light branch and on any unchanged original class branch. In a heavy branch the verifier constructs the full constraint-matrix transpose and checks y≥0, Aᵀy≥g, and the rational dual objective. It independently verifies the expected cutoff set and all382 records, so no endpoint or branch is omitted. The maximum class budget is kappa.

## Poisson summation and density

The zero mod2 class has no nonzero vector of length<2, so every such term has h<0. All2^18−1 other classes contribute at most2kappa, where the2 counts the two signs in each antipodal pair. Therefore

    rho (2pi/tau)^9 q(0) ≤ sum_(v in Lambda) h(v)
                          ≤ p(0)+2(2^18−1)kappa.

Multiplication by the packing-ball volume pi^9/(2^18 9!) gives the exact rational bound

    Delta18 ≤ (tau/8)^9 [p(0)+2(2^18−1)kappa]/[9!q(0)]
             <0.010343803448.

The alternative rho≤32768 is separately below this bound, proved using the rational pi upper bound22/7. Thus the conditional trace argument covers every lattice.
