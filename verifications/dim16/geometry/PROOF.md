# Full angular geometry of the frozen 16D log-hinge kernel

This proves only the angular kernel conditions. It is not by itself a
packing-density certificate or a spherical-count certificate.

The normalized angular polynomial Q is reconstructed exactly in dimension
16, with nonnegative Schoenberg coefficients and constant coefficient one.
Exact Sturm tests establish Q<=-1/20000 on [-1,1/2] and exactly one zero
in (1/2,1), with Q(1)>0. Thus Q has one transition from negative to
positive on the relevant angular interval. The bound |Q(t)|<=Q(1)=B on
[-1,1] follows from the normalized Gegenbauer addition theorem.

For v>=1, define phi(v) using the exact printed decimal coefficients in
fit.json: a nonnegative sum of v-1, (v-k)_+, log(v), and log(v/k)_+.
The linear knots are 1.02,1.05,1.1,1.2,1.4,2; the logarithmic knots are
1.2,1.4,1.5,1.6,1.65,1.75. Let h(v)=exp(-phi(v)). We have phi(1)=0.

## All-scale extension

The function g(v)=v phi'(v) is nondecreasing wherever differentiable,
with only nonnegative jumps at the knots. Indeed, a linear hinge
contributes v after its knot, and a logarithmic hinge contributes a
constant after its knot. Consequently, for z>=y>=1, both
h(lambda z)/h(lambda y) and h(lambda y)/h(lambda) are nonincreasing
in lambda>=1.

For a sorted triangle with side lengths 1<=u<=w, write
z=w^2=1+u^2-2ut and y=u^2. Its minimum-side-normalized cyclic expression is

    Q(t) + exp(phi(y)-phi(z)) Q((1-ut)/sqrt(z))
         + exp(-phi(z)) Q((u-t)/sqrt(z)).

The first angular value is negative. The other two angles are ordered.
If both remaining Q values are nonnegative, decreasing the two ratios
preserves nonpositivity. If the middle Q value is negative and the last
positive, first factor out the middle ratio. The bracket decreases with
scale because the quotient of the last and middle ratios decreases. While
the bracket is positive its product also decreases; after it becomes
negative, the total expression is already negative. If both Q values are
nonpositive, the claim is immediate. This proves propagation from minimum
squared radius one to every larger scale without differentiating at knots.

## Compact domain

The verifier starts with the exact rectangle
[-1,1/2] times [1,4] in (t,u). It discards only boxes on which 2tu>1
throughout, so they contain no sorted triangle. Every other box is
certified directly by interval arithmetic or by a centered mean-value
bound using interval derivatives. Accepted and discarded leaves have total
exact area 9/2, and are obtained solely by binary partition of the starting
rectangle. No floating acceptance test is used.

At a hinge the derivative is not unique. The implementation encloses both
one-sided derivatives. The profile and all composed expressions are locally
Lipschitz on this compact domain; the coordinatewise fundamental theorem
of calculus, using the almost-everywhere derivative, validates the same
mean-value bound. Points at a kink do not require an arbitrarily selected
classical derivative.

The lower bounds u^2>=1 and 1+u(u-1)+u(1-2t)>=1 hold on the entire
starting rectangle, even outside the sorted portion. They justify the
explicit interval intersections before square roots and logarithms.
The angle arguments are actual cosines for every such (t,u), so clipping
their polynomial range evaluation to [-1,1] is valid.

The completed replay accepts 31138 rectangles and discards 327 sorted-domain
exclusions, after 62929 tested boxes; the maximum binary depth is 46.

## Infinite tail

The exact final linear slope gamma is at least 7. For z>=y>=16,
phi(z)-phi(y)>=gamma(z-y) and phi(z)>=gamma(z-2).
If u>=4 and t>=-1/4, sortedness gives 1-ut>=1/2 and the middle cosine is
at most 1/u-t<=1/2. Only the last term can be positive, bounded by
B exp(-14 gamma)<1/20000.

If t<=-1/4, z-y=1-2ut>=3. The sum of the last two positive parts is at
most B[exp(-3 gamma)+exp(-14 gamma)]<1/20000. Exact Arb bounds verify
both comparisons. Thus the entire infinite tail is strictly negative.

## Replay and limitations

Run `python verify_geometry.py` with python-flint installed. It first runs
the exact structure/tail check, then the full compact covering. Assertions
are mandatory; both entrypoints reject Python -O and -OO. Precision is
256 bits. Nonfinite values, unresolved depth, elapsed budget, or leftover
boxes cause failure, never success. The saved partition is audit evidence;
replay recomputes all inequalities.

The radial Fourier function, its global signs, the final density pair
inequality, and any added spherical-count certificate require separate
proofs. They are not asserted here.
