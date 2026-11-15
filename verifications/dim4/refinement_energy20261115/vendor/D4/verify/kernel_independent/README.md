# Independent reconstruction of the V7 kernel

From the D4 folder:

```sh
python3 verify/kernel_independent/run.py --out /tmp/d4-independent
python3 verify/kernel_independent/negative_controls.py --out /tmp/d4-independent-controls
```

Both output directories must be new. The first command takes about 10 minutes
with one worker on an Apple M2 Pro; the recovered earlier run used 220 MB
peak RSS. `--workers 2` is supported. The controls take several minutes.
Dependencies are those of `../requirements.txt`, including a C++ compiler
for stage 1's fresh pair checks. No file outside the D4 package is an input.

`run.py` returns exit 0 only with `RESULT.json` containing
`stage: kernel_independent`, `status: PASS`, all **84 full and 28 exterior**
coefficient identities, all **7 pair pieces**, the exact origin/objective
identities, the component contracts and internal consistency checks.
Omitted cases, errors, wrong domains and mismatched coefficients fail.
`--only` is a diagnostic option and always fails the full coverage gate.

## What is independently computed

`kernel.py` reads the four `data/v7/*.json` files. It imports **no stage 1
or stage 5 code**. It uses FLINT integer matrix multiplication to assemble
the Gram blocks, its own Bernstein expansion, a three-term Chebyshev
recurrence, and an independent homogenization of the three rooted terms.
It implements each cap, tent, allocation, count, spherical, graph and
gauge component pointwise, including the origin constant.

The comparison side copies the stage 1 record into the output directory,
then runs its unchanged coordinator through the fresh scalar, component,
origin and seven pair proof gates. A worker adapter captures the exact
polynomial that stage 1 supplies to its cyclic subdivision routine. The
adapter skips that subdivision and returns UNKNOWN to the coordinator;
this process is an identity comparison, not a substitute sign certificate.
The coordinator's final coverage assertion therefore fails by design.
The driver independently requires all 112 captured comparisons to pass.

For the fixed package the 112 identities compare **1,059,688 nonzero
coefficients**. Polynomial equality is exact over the rationals, rather
than a numerical or finite-point comparison. The domains are compared too.

## Precise proof boundary

| Item | Evidence from this checker |
|---|---|
| Gram/count/signed-spherical cyclic parts | Exact coefficients on all 112 branch polynomials |
| Pair polynomial `D` and scalar polynomial `F(alpha r²)` | Exact coefficients on all 7 radial pieces |
| Origin `A`, raw objective and census objective | Exact rational equality; 49 origin equations |
| Actual stage 1 component parameters | 775 exact equations, including 49 radial basis polynomials |
| Nonpolynomial correction formulas and cap self term `J` | Written identity proof in [COMPONENT_IDENTITIES.md](COMPONENT_IDENTITIES.md), with exact parameter matching |
| Boundary diagnostics | 10 radial ownership, 16 allocation and 1371 graph checks |
| Nonpolynomial cyclic evaluations | 3728 sample points, **diagnostic only** |
| Kernel polynomial versus its own point definition | 238 sample comparisons, **diagnostic only** |

The written component proof handles irrational cap-diagonal roots and the
strict/closed allocation cutoffs. It does not infer identities from samples.
The interpretation of stage 1's interval routines remains a code-review and
written-proof obligation. This is not a second independent proof of all
continuum signs and is not an unconditional Lean verification of the concrete
kernel. The separate `../kernel_link/` check connects the stage 5 box formula
to stage 1 where the far-pair increment uses it.

The implementation is independent arithmetic and assembly code, but its
author inspected stage 1 to resolve conventions. Those choices are recorded
in [CONVENTIONS.md](CONVENTIONS.md); this is not a clean-room implementation.

## Failure controls and outputs

`negative_controls.py` requires:

* An unchanged three-case subset has equal pair/branch objects but fails
  full coverage, as intended.
* Perturbing a spherical weight in the independent implementation produces
  an exact mismatch.
* Changing one Gram-ray weight only in the independent input produces an
  exact mismatch.
* Ten malformed data variants are rejected, including nonintegral degree,
  denominator, cap count and gauge count; a fractional ray coordinate;
  a short origin row; a negative cap weight; duplicate count radii;
  a boolean dimension; and a floating scalar coefficient.
* Four malformed captured polynomial inputs are rejected before any integer
  conversion: a fractional coefficient, fractional or zero denominator,
  and a tensor with the wrong number of axes.

Its `CONTROLS.json` must have `status: PASS` and both perturbations marked
`FAILS_AS_REQUIRED`. These are tests of failure detection, not bound proofs.
`bonus_stage5.py` is an optional old pointwise diagnostic and is not part
of the exact identity claim.

`MANIFEST.json` identifies the local implementation and documents used by a
run. Stage 1's own manifest separately checks the copied record and the
original data. Manifest equality establishes artifact identity, not
mathematical correctness. The output contains `cases/` (one report per
branch, plus `pair.json`, `contracts.json` and `coordinator.json`), a complete
extraction log and the copied scratch record. No generated source or report
is written back into the publication package.
