# Exact link between the stage 1 and stage 5 kernels

Run from `preprints/D4`:

```
python3 verify/kernel_link/run.py --out /tmp/d4-kernel-link
```

The output directory must be new. Exit 0 and `RESULT.json` with `stage: kernel_link`
and `status: PASS` establish the identities below. This uses the same Python,
FLINT, NumPy and C++ dependencies as the other stages. It takes about two minutes
on an Apple M2 Pro, with one worker and under 0.3 GB RSS. Temporary working copies
are removed after success; `--keep-work` retains them for debugging.

## Statement

On both kernel-box domains

```
B_j = [1,R_j]^2 x [-1,-133/250],   j = 1,5,
R_1 = 501473/500000, R_5 = 102597/100000,
```

stage 5's explicit mathematical formula is the negative of the kernel assembled
and sign-checked in stage 1. This is an exact expression identity on the whole
domain, not a comparison of values at sample points. The stage 5 interval
program evaluates a lower envelope of that formula when kinks are undecided;
it is not claimed to evaluate the exact kernel on every interval.

There are 25 ordered radial cells. The 15 cells with `r`'s region at most `s`'s
region correspond to stage 1's exterior branches; exchanging the variables
covers the other ten. The first cell is `B_1`'s radial domain, so the `B_5`
comparison includes `B_1`. In this domain

```
|y-z|^2 = r^2+s^2-2rst >= 2+2(133/250) = 383/125 > 49/16.
```

Thus the third side is strictly outside the support `7/4`, including every
boundary of the box. Only the exterior branches are relevant here.

## What runs

1. `src/extract_stage1.py` copies the stage 1 source tree and data to a new
   working directory after checking its manifest. It runs the real coordinator's
   pair/component preflight. It captures the objects passed to the exterior
   branch checker, then calls the actual `build_signed` polynomial assembly for
   each of the 15 branches. The final triangle inequality subdivision is not run
   here; the full stage 1 run remains compulsory.
2. `src/stage1_terms.py` reconstructs the exact correction terms from those
   captured objects. It directly calls the original nested functions for the
   Bernstein features, cap threshold, tent argument and allocation distance.
   The remaining short formulas are stated below and tied to their source
   functions in that module. The original bound routines are also evaluated on
   375 singleton component cases. Those are diagnostics, not the identity proof.
3. `src/trace_stage5.py` executes the actual `kcert.minusK_lb` and
   `kformula.pieces` code with an exact symbolic number type. The interval
   primitives are replaced by their exact mathematical functions: positive part,
   clipping, closed indicators and the two named allocation capacities. The
   rational functions, exponentials and square roots in their arguments are
   preserved. Both traced expressions must agree exactly in each cell.
4. `src/compare.py` subtracts the stage 1 correction from this formula. Every
   nonpolynomial atom must cancel. The remaining polynomial is transformed to
   stage 1 coordinates and compared coefficient by coefficient with the freshly
   assembled stage 1 numerator. All 25 angular endpoint faces `t=-1` are checked
   separately. Altered constants, polynomial coefficients, a missing cell and
   an altered endpoint and a nonintegral denominator must all be rejected.

`CELLS.json` lists the 25 checked cells and exact coefficient counts.
`NEGATIVE_CONTROLS.json` records the rejected corruptions. Logs retain the fresh
stage 1 preflight and stage 5 trace. Serialized intermediate expressions are
created by the subprocesses during this run, never accepted from a supplied
cache, and deleted after success.

## Algebra of the comparison

Stage 1 uses

```
t = -1 + z ((r+s)^2-49/16)/(2rs).
K = C(r,s,z)/(dd (2rs)^18) - corrections.
```

`C` is an integer polynomial and `dd` a positive integer. The polynomial
comparison multiplies by `(2rs)^18`, which is strictly positive since `r,s>=1`.
It compares all coefficients in `Q[r,s,z]`. The atom-free difference must have
no denominator and angular degree at most 18 before this transformation.

The correction terms are:

- Caps: `weight (w(s) P(r,t)+w(r) P(s,t))/2`, with
  `c(r)=(r^2+d^2-a2)/(2dr)`. If `sigma=-1`, `P=[c+t<=0]`; if `sigma=1`,
  `P=[c-t<=0]`; otherwise, for `-1<t<1`,
  `P=clip((1-(c-sigma*t)/sqrt((1-sigma^2)(1-t^2)))/2,0,1)`.
- Tents: `weight (w(s) relu(1-2(r^2+d^2+2drt)) +
  w(r) relu(1-2(s^2+d^2+2dst)))/2`.
- Allocations: `weight (u(r)+u(s)) capacity(r^2+s^2+2rst)/2`, with
  `u(r)=(2r-2)^4` on these cells.
- Graphs: weight times the prescribed closed radial and angular indicators.
- Gauge: `(w(r) psi(s)+w(s) psi(r))/2`, where
  `w=-q/(1-psi)`, `psi` is the captured shell flag, and `q` is the complete
  unlifted pair expression extracted from stage 1, including the scalar
  exponential and the exact self-cap reductions.

The two-level capacity is 1 below `(8453/20000)^2`, 1/2 from that threshold
through `(11547/20000)^2` inclusive, and 0 above. The classical capacity is 0 at
and above 1/2, `1-2x` from 1/4 to 1/2, and
`min(1,2(1-sqrt(x))^2)` below 1/4. These are the functions used by both stages.

`src/symalg.py` represents a sum of canonical rational-function coefficients
multiplying products of function atoms. Equal atoms have the same function name,
parameters and exact argument. No algebraic independence assumption is needed:
identical expressions denote the same function. Simplifications using Arb are
performed only when its enclosure proves a branch decision on the entire cell.
Unresolved kinks remain atoms and must cancel structurally. An indicator's
argument can be scaled only by a factor proved positive on the cell.

## Boundaries and lower envelopes

Radial cells use `[1,R_1]`, then `(R_i,R_(i+1)]`. Count indicators are closed.
Every internal radial face belongs to the preceding cell; no continuity across
radial discontinuities is assumed. Stage 5's `dec` offset selects the branch
only: the polynomial is still bounded on the entire closed interval. The next
cell's omitted lower face is covered by the preceding cell.

For every smooth-cap trace, the denominator square is checked to be exactly
`c(1-t^2)` with rational `c>0`. Hence it is positive throughout `-1<t<=-133/250`.
The temporary rational interval used to select this branch does not restrict the
identity's domain; the polynomial identity proves positivity on the full open
angular interval.

At `t=-1`, the closed cap event for `|sigma|<1` is `P=[c(r)+sigma<=0]`.
The exact formula uses this event. The real stage 5 lower-bound routine uses 0
when its denominator square is not positive. A separate symbolic comparison
checks that its endpoint expression is the exact formula minus precisely those
caps. Their coefficients are nonnegative (checked weights and Bernstein
coefficients on the radial basis interval), so this omission is a valid lower
bound, including a zero cap threshold. The singleton face at `r=1` and every
upper radial face is already included by the same symbolic identities.

For intervals crossing a kink, the existing stage 5 proof of lower-envelope
soundness remains required: [kernel_box_method.md](../stage5_kernelbox/proofs/kernel_box_method.md),
section 3. In brief, positive part and cap probabilities admit zero as a lower
bound, clipping is at least `max(0,inf x)` when below 1, and both allocation
capacities are nonincreasing. Every such term has a nonnegative multiplier in
`-K`; self terms in the gauge are decided and asserted so by the real code.
This symbolic identity check does not reprove Taylor's theorem or the Arb
implementation.

## Trust and scope

The proof combines exact executable algebra with the written meanings of the
small helper functions above. The symbolic adapters for nonpolynomial terms are
reviewable translations, not a formal verification of Python source semantics.
This is a link between two existing implementations, not a wholly independent
implementation of every component; see [kernel_independent](../kernel_independent/)
for the separate assembler. Python, FLINT exact arithmetic and Arb enclosures
remain in the computational trust base. Neither the singleton diagnostics nor
source hashes establish mathematical equality. The equality gate is the exact
coefficient/atom comparison, and the interval sign claims still require the full
stage 1 and stage 5 verification runs.
