# Dimension-parametric profile optimizer and exact certificates

## Scope and status

`pilot.py` is a numerical discovery tool, not a proof. `certify.py` produces rational interval and profile certificates; `verify_certificate.py` replays them using Python and SymPy, without a numerical optimizer. The verifier is conditional on the stated arithmetic trace-gap theorem for the chosen dimension. That theorem and its prior density assumption must be supplied separately. The low-density branch of that theorem must be combined with the certified high-density branch by taking their maximum.

This package publishes verified lattice-only upper bounds, without an exhaustive worldwide-priority claim. Per-dimension `prior_certificate.json` files for 14D, 15D and 18D were fetched from AstroNinja389/lattice-packing-bound, certificates/nNNN.json, and their claims are exactly replayed. Source-file bounds are named comparison baselines; see the separate lattice-results table for literature qualifications.

## Inputs that really depend on dimension

* Dimension n, Gaussian scale tau, squared-radius partition E, polynomial degree K.
* A proved implication: if rho exceeds a specified threshold, any n distinct opposite pairs in one nonzero mod-2 class have sum of squared radii at least 2n+delta.
* The prior upper bound used to prove that implication, and its exact verification.

The polynomial discovery code uses alpha=n/2-1, not a fixed integer alpha. It uses gamma(n/2+1) for exploratory density reporting. Exact verification uses rational half-integer binomials for odd n. Gaussian scales differ in the provided sources (e.g. 12D tau14,14D tau16,18D tau20); the cutoff for the polynomial tail is always tau*E[-1].

## Exhaustive first-n profile disjunction

Long-bin variables count opposite pairs, not individual vectors. Counts are relaxed to nonnegative reals. A class has either no short pair, or exactly one pair of squared norm below2. For a short bin j, its squared norm is bounded above by E[j+1]. All long-bin radius excesses also use upper endpoints, never lower endpoints.

For an empty class set count=n and short=0; for a short class set count=n-1 and short=E[j+1]-2. Profiles with fewer than count long pairs use the greedy monotone bound truncated at count-1. Otherwise enumerate h, the bin containing the count-th long pair. The necessary constraints are:

* cumulative counts bounded by the ordinary relative/anchored spherical-code caps;
* at most count pairs strictly before h;
* at least count pairs at or before h;
* sum_(b<h) (excess[h]-excess[b])*x[b] <= count*excess[h]+short-delta.

These are necessary constraints on the first n pairs, not the average radius of all pairs. Adding very distant negligible-weight pairs therefore cannot evade the trace-gap condition. Every feasible cutoff is included. A rational nonnegative LP dual bounds each cutoff. The verifier checks all dual inequalities and exact coverage of the cutoff cases, including the empty class. No dimension-specific weighted-empty-class inequality is needed.

## Exact function and objective

p(x)=sum d_k L_k^(n/2-1)(x)/L_k^(n/2-1)(0), q(x) has coefficients (-1)^k d_k in the same basis. The physical function is h(v)=p(tau*|v|^2)*exp(-tau*|v|^2/2). The verifier reconstructs p and q from d, proves q>=0 on [0,infinity) and p<=0 beyond tau*E[-1] by exact root counts, and verifies rational Bernstein interval majorants for every bin. A positive truncated exponential series supplies rational upper bounds for all bin heights. Long heights are replaced by backward maxima.

For R=(p(0)+2*(2^n-1)*kappa)/q(0), the density upper expression is (tau/8)^(n/2)*R/Gamma(n/2+1). For even n it is rational. For odd n=2m+1, its square times pi is the rational number

(tau/8)^n * [4^(m+1)*(m+1)!/(2m+2)!]^2 * R^2.

A rational reported upper U is certified by U^2*pi_lower being at least this number. The lower bound for pi is obtained from Machin's identity using a four-term lower alternating sum for atan(1/5) and the upper bound atan(1/239)<1/239. Thus neither floating gamma nor floating pi enters a final certificate.

## Numerical controls and failure handling

K61 pilots use HiGHS primal/dual feasibility tolerances1e-9, normalized Gaussian-Laguerre bases, and positive normalized1e-8 Fourier/tail margins. Initial grids alone are insufficient: candidates with negative extrema are rejected. `check_sign.py` searches continuous extrema with double precision and reevaluates them at70 digits to seed denser neighborhoods. This is only a diagnostic. Exact global root checks decide acceptance. Fourier positivity failed for an early18D candidate even with positive grid margins; the candidate was discarded and regenerated.

`REFINE=8` subdivides each source bin into eight bins. Always compare a trace-gap pilot against a no-gap pilot on the same refined partition; otherwise ordinary partition improvement is incorrectly credited to the arithmetic theorem.
