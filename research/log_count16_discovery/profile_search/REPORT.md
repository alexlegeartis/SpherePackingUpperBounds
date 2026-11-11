# Fixed-Q nested log-hinge experiment, 16D

**Status: bounded discovery experiment completed; no new packing bound.**

The single new fit and single new radial LP give sampled expression
**0.0249281075589411**, versus the frozen old-cone dual refit
**0.0249333317501032** and matched quadratic control **0.0249421111607885**.
The new expression remains **0.00000689755894108 above 0.02492121**.
It fails the requested promotion threshold, so no full proof or second fit was
attempted. This folder must not be cited as a density certificate.

## Scope and reproduction

- Fixed exact angular input SHA-256:
  `e6f6d85ca4fb160536f3ebc74d94f9ab7ea659793cad94258c8a089447a3f40b`.
- The copied `experiment.py` and `inputs/angular_dim16.json` are unchanged.
- `frozen_control/` contains immutable copies of the source profiles, audits,
  exact-Q/tail result, radial controls, radial rows, and source-profile-0 duals.
- One SLSQP fit, 57 iterations, 51.44 seconds through fitting. Total fitting,
  separation, and radial solve: **63.45 seconds**.
- One degree-63 candidate LP, scale `s=2`, tail cut `X=40`, 13,606 rows.
  Cached controls use precisely the same pair/Fourier/sign sample locations,
  normalization, objective, row-scaling rule, tolerance, and solver options.
  The controls were not re-solved. Their own global sign audits had failed.
- No remote repository writes, no denied actions, and no access blockers.

Run from the workspace:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/density_log_hinge_profile16_20261111/log_hinge.py
OPENBLAS_NUM_THREADS=1 python research/density_log_hinge_profile16_20261111/summarize.py
```

This reproduces the bounded numerical experiment, not a global certificate.
Main artifacts: `fit.json`, `cyclic_audit.json`, `radial_result.json`,
`summary.json`, `run_status.json`, `candidate_pair_dual_pricing.npz`, and
`manifest.json`. Source and result hashes exclude the manifest itself and the
live redirected summary log.

## Nested profile cone

The original seven nonnegative features are preserved exactly, in order:

`v-1`, followed by `(v-tau)_+` at `tau = 1.02, 1.05, 1.1, 1.2, 1.4, 2`.

Seven nonnegative features are appended:

`log(v)`, followed by `(log(v/sigma))_+` at
`sigma = 1.2, 1.4, 1.5, 1.6, 1.65, 1.75`.

The first seven feature arrays are byte-for-byte identical to the old basis.
All-zero added coefficients recover the control symbolically. Expanded-width
floating-point matrix multiplication changes evaluated phi by at most
5.68e-14 in the nesting check; this is summation roundoff, not a feature change.

The total coefficient of the seven linear features is constrained to be at
least 7. This includes the old control and retains a Gaussian tail even when
the optimizer shifts all earlier decay into logarithmic features.

Ignoring coefficients below 1e-13 only for this human-readable display, the
new candidate is approximately:

```
phi(v) = 7 (v-2)_+
       + 11.9605832710351 log(v)
       + 6.69847004072951 (log(v/1.2))_+
       + 2.67484278860294 (log(v/1.4))_+
       + 16.1696638052215 (log(v/1.5))_+
       + 39.1289146621932 (log(v/1.6))_+.
```

The actual saved coefficients, not this rounded display, were used by the LP.
Their exact decimal sum over linear features is greater than 7 and all are
nonnegative. Log hinges at 1.65 and 1.75 are zero in this fit.

## Pricing and separation

The fit minimizes the actual tangent objective recovered from the source
profile-0 radial LP. The dual array was not mislabeled as belonging to the
later old-cone refit. Its raw price for a profile perturbation is

`(-scaled_marginal / fixed_row_scale) * 2*C*exp(X)/B * h*(1-h)`.

That price is contracted with the full new 14-feature vector. The feasible
starting point is the better frozen old-cone refit with all seven added
coefficients zero. The normalized tangent objective fell from
7.49548586088725 to 7.15603368362682, about 4.53%.

All constraints retain all three **signed** Q contributions. There is no
positive-part surrogate. Coverage:

- 343,682 pooled training rows, including every old separated witness;
- 299,542 dense audit rows;
- 240,000 independent seeded random audit rows;
- both `u^2=tau` and `w^2=tau` fibers at every old and new knot;
- a finer mesh around the sharp isosceles obstacle;
- 150 bounded local separation runs.

Worst separated cyclic value: **-1.99717962178e-7**, at approximately
`(t,u)=(0.181044622802,1)`. This is numerical feasibility, not a compact-domain
proof. The analytic parameter Jacobian agrees with finite differences to
1.58e-11 in the recorded spot check.

## What bottleneck changed

At the old primary priced radius `v=1.61532400565`, phi's slack above the
isosceles lower obstacle fell from **0.629190** to **0.064330**.
The new radial dual shifts to:

- `v=1.63678879654`, tangent weight 0.0147815, obstacle slack **0.000768739**;
- `v=1.64036626169`, tangent weight 0.0133628, obstacle slack **0.00431359**.

There is also a much smaller price near `v=1.94087333422`; see `summary.json`.
The fit spends more decay after the obstacle: at the old secondary radius
`v=1.72980289042`, its obstacle slack increases from 4.01217 to 5.64547.

The new cone therefore removes much of the previously wasted decay, but the
dominant current prices now sit close to the same fixed-Q isosceles obstacle.
This is a diagnostic, not an optimality or impossibility proof.

## Structural rationale, conditional on the compact cyclic gate

For nonnegative coefficients in this cone, `g(v)=v*phi'(v)` is nondecreasing:
between knots it has nonnegative derivative equal to the active linear slope;
each linear hinge adds an upward jump `tau*a`, and each log hinge an upward
jump equal to its coefficient. Thus ordinary convexity of phi is unnecessary.

For a normalized sorted triangle with squared sides `1 <= y <= z`, define

`A(lambda)=phi(lambda*z)-phi(lambda*y)` and
`D(lambda)=phi(lambda*y)-phi(lambda)`.

Both are nondecreasing for `lambda>=1`, since their derivatives multiplied by
lambda are differences of the nondecreasing function g. The signed cyclic
expression is

`q0 + exp(-A(lambda)) * [q1 + exp(-D(lambda))*q2]`.

The frozen Q has negative cap and one sign transition. Ordered triangle angles
therefore imply that `q2<=0` forces `q1<=0`, when the expression is immediately
negative. If `q2>0`, the bracket is nonincreasing. Once nonpositive it remains
safe; while positive, multiplication by the nonincreasing exponential cannot
increase it. Consequently a **proved** normalized compact cyclic inequality
would extend to all scales. That compact proof has not been supplied here.

The linear-tail floor also avoids a pure-log asymptotic failure. For `v>=2`,
`phi(v)>=7*(v-2)` and `phi'(v)>=7`. For `u>=4`:

- If `t>=-1/4`, the middle angle lies in Q's negative cap, and the only
  potentially positive contribution is at most `B*exp(-98)`.
- If `t<=-1/4`, then `z-y>=3`, so the two potentially positive contributions
  sum to at most `B*(exp(-21)+exp(-98))`, about 6.32e-6, below the frozen
  cap margin 1/20000.

This elementary structural tail argument does not certify the radial Fourier
or physical-tail polynomials. Those global radial sign gates, the annular
pair inequality, and the compact cyclic interval gate were **not run** because
the sampled expression already fails the requested comparator threshold.
