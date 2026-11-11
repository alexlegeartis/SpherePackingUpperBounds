# Discovery package verification, 2026-11-11

All tests below passed with Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0,
CVXPY 1.9.3, Clarabel 0.11.1, and python-flint 0.9.0. No dependency
installation, fresh optimization, remote write, or integration edit was
performed as part of this packaging task.

## Deterministic and rejection tests

`python test_discovery.py`: **7 tests passed**.

- All eight Python files compiled. The six CLI entrypoints returned help.
- Every Python file rejected optimized execution with `python -O` before
  proceeding with its work.
- The final saved numerical checkpoint reconstructed the exact source
  byte-for-byte, SHA-256
  `235920216afbd10b38be5efd3672e9d94b5ec5af8fc54aca3a63d77d9997721f`.
  All 13 off-origin shifts were zero.
- Existing output paths were rejected without changing their contents.
- The search dry-run contained all required final-stage settings, defaulted
  to two rounds and `wedge5.json`, and resolved the budget and far boxes to
  bundled local paths. Zero rounds were rejected. No optimizer was launched.
- A matrix requiring an off-origin shift was rejected without writing a
  prepared output.
- Copies of both packages relocated to a temporary directory, invoked from
  a separate working directory, produced canonical exports identical to the
  sibling proof package for all three payloads: witness, exact PSD witness,
  and prerequisites. Export-only correctly reported incomplete verification.
- Mutated budget bytes were rejected before creating output in both the
  exporter and search wrapper. The uncertified numerical checkpoint was
  rejected as input to the frozen complete-replay command.

Additional input checks confirmed byte identity with the original warm start,
final checkpoint, exact source, cyclic fiber rows, far boxes, and old centroid
budget. The dimension-five budget has caps 83/8/5/3/2/1 and 104 bins. No newer
spherical budget was substituted.

## Complete final proof-suite replay

Ran `python export_and_check.py --output <new temporary directory>` against
the final sibling `lens_verification_package`. It re-exported all inputs,
confirmed canonical equality, copied the final proof scripts, and reran
complete verification in the fresh output directory.

The sibling proof manifest SHA-256 was:

```text
d76a9cacc91a8acc9c481b904d6289e18a64cdf7635f3f905a9f8c9498ec3840
```

The replay finished with `FROZEN_EXPORT_AND_COMPLETE_VERIFICATION_PASS` and
`complete_verification_passed: true`. It checked the final hardened algebra,
budget query, and far-verifier scripts, not the earlier provisional copies.

- Frozen source, exact payload, centroid geometry, and boundary conventions:
  passed.
- Exact rational positive definiteness and Schur factorization: passed.
- Radial/Fourier global signs, zero Gaussian corrections, and outward objective:
  passed.
- Exact independent algebra checks, all step budgets, candidate binding, and
  explicit half-open count endpoints: passed.
- Pair verifier: 121 checked boxes, maximum depth 11, no failure.
- Cyclic verifier: 115 checked boxes, 52 passed, 6 excluded, maximum depth 14,
  no failure.
- Far verifier: 2,879 checked boxes, 915 passed, 526 excluded, 11 relaxed boxes,
  maximum depth 22, no failure.

The final outward density enclosure was:

```text
[0.5099961464057082278011891612604254673934488294265848385588116530194451139022 +/- 3.91e-77]
```

The verified strict upper is **0.509996147**. Scope: unrestricted sphere
packings, with the published near-neighbor count bounds documented by the
sibling proof package as external mathematical inputs.

## Limits

The test run did not rerun the numerical optimizer or reconstruct the complete
exploratory search history. The saved final bounded warm-start stage, its
inputs, and command are supplied for that purpose. Numerical outputs may vary
and remain uncertified until separately prepared, completely verified, and
reviewed. The complete replay here applies only to the frozen exact source.
