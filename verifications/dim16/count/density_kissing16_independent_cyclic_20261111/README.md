# Exact widened spherical count certificate in dimension 16

**Result:** the frozen witness certifies that a spherical code in S^15 with every distinct inner product at most **5003/10000** has at most **7802** points.

The real upper bound has the exact rational enclosure

7802.215167128723216019597445 < U < 7802.215167128723216019597446.

This is a spherical counting certificate. A density theorem additionally requires its own geometric/profile/radial hypotheses, which this directory does not verify.

## Frozen inputs

Both are read-only inputs in `../density_kissing16_reweight_frozen_20261111/`:

- `candidate.json`: SHA256 `b3225ad4cce819395ce3bd47613f15697af4646fb722ded1916020b5e52ba573`.
- `factors.json`: SHA256 `6860eca2c82fdb1b1ebb5d3a4108a216d46c900d6fac7bdeb0d880272e027954`.

All printed decimal numbers, including candidate numbers, are interpreted as exact rational decimals. The source factor orientation is X=L^T L+lambda I. Nothing in the frozen inputs was modified. The old source contact SOS identity is not assumed or replayed: the new pair and cyclic inequalities are checked directly.

## Independently reconstructed witness

Use the normalized Gegenbauer polynomials P_j^m(1)=1 with

(k+m-2) P_(k+1)^m(x) = (2k+m-2)x P_k^m(x) - k P_(k-1)^m(x).

The homogeneous transverse factor Q_k is the polynomial obtained from P_k^15 with g=t-uv and z=(1-u^2)(1-v^2): Q_0=1, Q_1=g, and

(k+13)Q_(k+1)=(2k+13)g Q_k-kz Q_(k-1).

Set p_k(u)=(P_0^(16+2k)(u),...,P_(16-k)^(16+2k)(u)) and

K(u,v,t)=B_11+sum_k theta_k Q_k(u,v,t) p_k(u)^T F_k p_k(v).

The extra nonorigin kernel is R(t)=sum_l r_l P_l^16(t), with zero origin and cross components. The global function is f(t)=1+sum_l a_l P_l^16(t).

All theta_k, a_l, r_l are nonnegative exact rationals; all source lambdas are positive. Thus every reweighted local factor remains PSD. The Q addition theorem gives positive semidefiniteness on the orthogonal projections in dimension 15, including degenerate u=+/-1 by the polynomial extension. The scalar features p_k(u) preserve PSD, as do B and the separate nonorigin R kernel. The origin value is

C=B_00+theta_0 1^T F_0 1.

## Pair gate: PASS

The separate directory `../density_kissing16_independent_pair_20261111/` contains the exact PSD, objective and pair proof, its complete rational Sturm chains, and a Python-integer/Fraction replay.

The degree-32 pair polynomial

f(t)+2[B_01+theta_0 1^T F_0 p_0(t)]+K(t,t,1)+R(1)

is strictly below -9/10000000 throughout the full closed interval [-1,5003/10000]. Both the original and margin-shifted Sturm chains have endpoint variation counts (17,17), hence no roots, with negative endpoint values. Two direct pair reconstructions agree exactly.

## Cyclic gate: PASS

The polynomial

H(u,v,t)=K(u,v,t)+K(u,t,v)+K(v,t,u)+R(u)+R(v)+R(t)

has 2673 nonzero rational monomials, total degree 32, and degree at most 16 in each variable. Its full permutation symmetry is checked exactly.

The domain is the full cube [-1,5003/10000]^3 intersected with

G(u,v,t)=1+2uvt-u^2-v^2-t^2 >= 0.

The Bernstein proof starts from exactly transformed rational degree-(16,16,16) coefficient tensors for both H and G. Dyadic boxes form a complete partition of the normalized unit cube. Every terminal box is justified by exactly one of:

1. All Bernstein coefficients of H are strictly negative.
2. All Bernstein coefficients of H+lambda G are strictly negative for a specified nonnegative rational lambda. Since G>=0, H<=H+lambda G<0.
3. All Bernstein coefficients of G are strictly negative, so the box is infeasible.
4. The box lies strictly outside u<=v<=t. Exact symmetry covers the discarded chambers. Boundary equality is never discarded by this test.

The partition has 2155 traversal events: 1077 splits, 641 direct-negative leaves, 259 Gram-multiplier leaves, 107 Gram-negative leaves and 71 symmetry exclusions. Maximum total depth is 46. The exact sum of all terminal box volumes is 1, with no unfinished boxes.

`bernstein_gram.py` discovers and verifies this partition using outward-rounded binary64 intervals. Its positive minimum acceptance margin is approximately 5.30e-9. This numerical margin is not needed by the final exact replay.

`replay_exact.py` reconstructs the rational Bernstein tensors and verifies **every saved proof event entirely by arbitrary-size integer arithmetic**, including every leaf coefficient sign. Initial H and G share a common denominator with 302 decimal digits. Each dyadic subdivision uses integer de Casteljau matrices and multiplies the common denominator by 2^16. Each multiplier is stored by its exact integer numerator and denominator. The exact replay passes in approximately 17 seconds. No floating-point/Arb evaluations or sampled points are used for proof acceptance in that replay.

The independent audit in the pair directory separately reconstructed all 2673 polynomial coefficients from the bilinear factor representation, checked all 4913 initial interval enclosures against its own rational transform, and checked 29478 first-subdivision interval coefficients. Its further independent integer-grid traversal also passes: a separately reconstructed rational tensor is rounded upward to a 2^-80 integer grid, then every midpoint, multiplier and sign test is verified by integer arithmetic. The exact partition hash matches. It certifies the uniform cyclic upper bound -1804896546625092345704439846467/2^128, approximately -5.30411423594e-9. This is conservative; it is not a claim that the true maximum is close to zero.

## Why these gates give the count

For a code of N points, positive definiteness of f with constant coefficient 1 gives sum_(i,j) f(<x_i,x_j>) >= N^2.

For each center i, sum the positive semidefinite local kernel over the origin and all other code points. Summing over centers gives a nonnegative total. The pair inequality bounds all its origin/cross/diagonal terms by NC minus the ordered off-diagonal f-sum. Every unordered triple contributes twice the full cyclic sum, which is nonpositive by the certified H inequality. Therefore

sum_(i!=j) f(<x_i,x_j>) <= NC,

and N^2 <= N[f(1)+C]. For N>0, N<=U=f(1)+C. Since N is an integer and U<7803, N<=7802.

## Portable top-level audit

Run `python density_kissing16_independent_cyclic_20261111/audit_count.py` from any directory, using a Python environment with NumPy and python-flint installed. Preserve these three sibling directories beneath any parent:

- `density_kissing16_reweight_frozen_20261111`
- `density_kissing16_independent_pair_20261111`
- `density_kissing16_independent_cyclic_20261111`

`DEPENDENCIES.json` specifies all source-directory dependencies and tested package versions. No discovery code, original SDP archive, network, or fixed workspace is required. SciPy and SymPy are not required for the complete audit.

The top-level audit rebuilds the pair and cyclic polynomial certificates, replays the Sturm chains, independently reconstructs the cyclic polynomial, runs both integer cyclic replays, verifies the original frozen input and mathematical proof-event hashes remain unchanged, tests explicit `-O` rejection for every Python entrypoint, tests top-level `-OO` rejection, and refreshes `CERTIFIED_COUNT.json`, `AUDIT_RESULT.json`, and manifests. Audit results are recomputed, not accepted from existing logs.

Every Python entrypoint now raises an explicit RuntimeError before any verification if assertions are disabled. Do not use `python -O` or `python -OO`. The full six-stage audit passed both in the original location and after relocation of only the three dependency directories beneath a parent not named `research`, with a different working directory. All 11 Python entrypoints reject `-O`, and top-level `-OO` also rejects; see `AUDIT_RESULT.json` and `PORTABILITY_RESULT.json`.

## Replay

From the workspace root, use the installed single-thread Python:

`PY=research/density_solid_cyclic_20261111/.venv/bin/python`

- `$PY research/density_kissing16_independent_pair_20261111/verify.py`
- `$PY research/density_kissing16_independent_pair_20261111/replay_chain.py`
- `$PY research/density_kissing16_independent_cyclic_20261111/bernstein_gram.py 200000`
- `$PY research/density_kissing16_independent_cyclic_20261111/replay_exact.py`
- `$PY research/density_kissing16_independent_pair_20261111/audit_cyclic.py`
- `$PY research/density_kissing16_independent_pair_20261111/replay_gram_integer.py`

The Gram proof always rebuilds exact tensors from the frozen inputs. The earlier `bernstein.py` run without Gram multipliers stopped at a 100000-box work limit and is retained as an explicitly incomplete experiment. Its optional cache reuse is not used by the final certificate.

## Numerical challenges, not proof

`challenge.py` tested 1642453 new random, boundary and special-locus triples, local SLSQP searches, and three differential-evolution searches. `challenge_edges.py` additionally tested 1333668 points on the exact Gram-boundary/cap edges and refined their numerical maxima. No positive point was found. `check_points_arb.py` checked selected rational and algebraic feasible points at 192-bit precision.

A formerly missed feasible boundary value is about -6.0357861e-7, larger than the original discovery separator's value but still strictly negative. A separate optimizer point had tiny negative Gram determinant and was explicitly rejected as a feasible proof point. Neither sampling nor any Gram-invalid optimization point is used for certificate acceptance.
