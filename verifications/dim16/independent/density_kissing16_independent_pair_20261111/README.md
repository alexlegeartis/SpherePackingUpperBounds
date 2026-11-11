# Independent exact pair verification, dimension 16

## Result

For the frozen candidate with cap `5003/10000`, the exact reconstructed pair
polynomial `H` has degree 32 and satisfies

`H(t) < -9/10000000` for every `-1 <= t <= 5003/10000`.

All candidate weights and Gegenbauer coefficients are exact nonnegative
rationals. All 34 source blocks have the form `L^T L + lambda I`, with strictly
positive exact `lambda`. Thus the fixed `B` block and all reweighted local
blocks are PSD. The source factors' exact original objective also reconstructs
correctly as an orientation/indexing cross-check.

The new objective, **not yet a count bound**, is enclosed rigorously by

`7802.215167128723216019597445 <= U <= 7802.215167128723216019597446`.

The exact rational value is retained in `result.json` and
`pair_coefficients.json`.

**The pair proof alone does not check the cyclic inequality. The additional
independent integer cyclic replay described below does check the full Gram
feasible domain. Applying the spherical theorem to infer a count or density
bound remains outside this folder's scope.**

## Frozen inputs

The verifier refuses inputs whose hashes differ from:

- `../density_kissing16_reweight_frozen_20261111/candidate.json`:
  `b3225ad4cce819395ce3bd47613f15697af4646fb722ded1916020b5e52ba573`
- `../density_kissing16_reweight_frozen_20261111/factors.json`:
  `6860eca2c82fdb1b1ebb5d3a4108a216d46c900d6fac7bdeb0d880272e027954`

JSON decimal coefficients, including exponent notation, are interpreted exactly
as printed. No binary floating-point value enters the mathematics. Runtime
measurements are the only floating-point output of this verifier.

## Independent reconstruction

`verify.py` does not import the discovery implementation or another verifier.
It reconstructs the formula from the factor encoding and the definitions used
in `../density_kissing16_joint_repair_20261111/{screen,reweight}.py`.

For normalized Gegenbauer polynomials in dimension `m`, it uses

- `P_0=1`, `P_1=t`;
- `P_(j+1)=((2j+m-2)t P_j - j P_(j-1))/(j+m-2)`.

It checks `P_j(1)=1` and `P_j(-1)=(-1)^j` exactly. Global `f` and nonorigin-only
`R` are reconstructed from all 24 stored coefficients:

- `f(t)=1+sum_(j=1)^24 a_j P_j^16(t)`;
- `R(t)=sum_(j=1)^24 r_j P_j^16(t)`.

For local index `k`, the vector `p_k` contains degrees `0,...,16-k` in dimension
`16+2k`. The frozen local block is `F_k=L_k^T L_k+lambda_k I`.
The homogeneous `Q` recurrence in dimension 15 is evaluated at `(t,t,1)` using
`g=1-t^2`, `z=g^2`, and

`Q_(k+1)=((2k+13)g Q_k - k z Q_(k-1))/(k+13)`.

Every `Q_k` is checked to equal `(1-t^2)^k` exactly. The pair polynomial is

`H(t)=f(t)+2 B01+B11+2 theta_0 1^T F_0 p_0(t)`
`     +sum_(k=0)^16 theta_k Q_k(t,t,1) p_k(t)^T F_k p_k(t)+R(1)`.

For every local block, the quadratic expression is computed independently in
factor form (squares of `L p`) and dense-matrix form; their exact polynomials
are compared. The cross expression is likewise checked in both forms.

The origin term and objective are

- `C=B00+theta_0 1^T F_0 1`;
- `U=f(1)+C`.

`R` applies only to the nonorigin block and contributes no origin term.

## Exact Sturm proof

For `H` and separately for `H+9/10000000`, an exact rational Sturm sequence is
computed. Each element is reduced to a primitive integer polynomial only by a
**positive** scalar, so signs are preserved. Both sequences have degrees
`32,31,...,0`, and both have 17 sign variations at each endpoint. Thus there
are no roots in the open interval. Endpoint values are strictly negative.
Continuity proves strict negativity on the whole closed interval. The shifted
proof yields the uniform slack stated above.

Both full primitive integer Sturm chains are retained, not merely their sign
vectors. `replay_chain.py` rechecks all chain identities and endpoint signs
using only Python integers and `Fraction`, independently of FLINT/SymPy. Its
pseudo-division always multiplies by the absolute value of the divisor's
leading coefficient and removes positive integer content, preserving the
remainder sign. It also verifies input and certificate hashes and the exact
objective enclosure.

## Reproduction

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
research/density_solid_cyclic_20261111/.venv/bin/python \
research/density_kissing16_independent_pair_20261111/verify.py

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
research/density_solid_cyclic_20261111/.venv/bin/python \
research/density_kissing16_independent_pair_20261111/replay_chain.py
```

To reproduce the additional independent cyclic checks, run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
research/density_solid_cyclic_20261111/.venv/bin/python \
research/density_kissing16_independent_pair_20261111/audit_cyclic.py

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
research/density_solid_cyclic_20261111/.venv/bin/python \
research/density_kissing16_independent_pair_20261111/replay_gram_integer.py
```

The first pair script uses exact `python-flint` rational arithmetic with
`ctx.threads=1`. The pair Sturm replay uses only the standard library.

## Artifacts

- `verify.py`: independent factor/polynomial reconstruction and Sturm generation.
- `pair_coefficients.json`: exact monomial coefficients and exact objective.
- `pair_sturm_chain.json`: full unshifted primitive integer Sturm chain.
- `pair_plus_9e_minus_7_sturm_chain.json`: full shifted chain.
- `result.json`: hashes, exact values, endpoint signs, and conclusions.
- `replay_chain.py`: separate standard-library certificate replay.
- `replay_result.json`: replay verification result.
- `run.log`, `replay.log`: execution output.

## Additional independent cyclic-construction audit

`audit_cyclic.py` independently expands the cyclic kernel directly from the
factor bilinear forms, rather than the other implementation's dense basis
matrix. All 2,673 rational coefficients match
`../density_kissing16_independent_cyclic_20261111/polynomial.json` exactly.
It then performs affine substitution and power-to-Bernstein conversion as
separate rational transformations and checks all 4,913 saved binary64 coefficient
intervals by exact rational comparison. Finally, every child coefficient from
the first split on each of the three axes is compared against exact rational
de Casteljau values (29,478 intervals total). All checks passed.

This audit is tied to the script/polynomial/cache hashes in
`cyclic_audit_result.json`. It is a construction and rounding check, not a
claim that the global cyclic covering search completed.

## Independent integer replay of the completed cyclic covering

`replay_gram_integer.py` is a separate implementation of the successful Gram
multiplier traversal. It uses the independently reconstructed exact tensor in
`exact_cyclic_bernstein.json`. Each rational initial coefficient is rounded
**upward to an integer multiple of `2^-80`**. For subdivision, exact integer
`ceil((a+b)/2)` is applied at every de Casteljau step. These are upper bounds
by induction, without floating-point assumptions.

Floating-point calculations reproduce only the original split choices and
proposed multipliers. A direct leaf is accepted only when the integer upper
coefficients of `H` are all negative. A Gram-excluded leaf requires all integer
upper coefficients of `G` to be negative. For a multiplier `lambda=p/q>=0`,
the integer coefficients `q*upper(H)+p*upper(G)` must all be negative. Thus
`H+lambda G<0`, and `H<0` wherever `G>=0`. The multiplier's exact dyadic
rational value is retained. Sorted-chamber exclusion uses exact rational
endpoints and strict inequalities.

The traversal must exhaust its stack, its exact box volumes must sum to one,
and its counts and ordered partition hash must match the other implementation.
Every terminal leaf is retained in `gram_replayed_leaves.json`. The result is
recorded in `gram_integer_replay_result.json`; the script asserts all checks
before writing a passing result. The fixed-point representation makes this a
separate exact proof arithmetic path rather than another floating-point audit.

The completed integer replay checked 2,155 nodes, with 641 direct and 259
multiplier leaves, 107 Gram exclusions, and 71 symmetry exclusions. The exact
covered volume is one, maximum depth is 46, and the partition hash matches
`9989c48afc1068a38a21c04597604cfa6b3927b015dcd364957c12c317769576`.
Its exact minimum cyclic margin exceeds `1/200000000`, proving the cyclic
polynomial is strictly below `-5e-9` throughout the Gram-feasible domain.

## Portable complete audit and optimization guards

Every Python audit entrypoint explicitly raises RuntimeError when `__debug__` is false. Neither `-O` nor `-OO` is permitted. The complete top-level count audit is `../density_kissing16_independent_cyclic_20261111/audit_count.py`; it rebuilds and replays all certificates, tests guard rejection, and refreshes manifests. Its `DEPENDENCIES.json` specifies the three sibling directories and Python dependencies. Full normal audits passed both in the original location and after relocating only those directories beneath a parent not named `research`. See `PORTABILITY_RESULT.json` in the cyclic directory.
