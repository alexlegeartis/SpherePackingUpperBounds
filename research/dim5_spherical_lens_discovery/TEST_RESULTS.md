# Spherical-lens discovery package verification, 2026-11-11

Tests used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, CVXPY 1.9.3,
Clarabel 0.11.1, and python-flint 0.9.0. No fresh optimization, dependency
installation, remote write, or integration edit was performed.

## Deterministic and rejection checks

`python test_discovery.py`: **7 tests passed**.

- All eight Python files compiled; all six CLI entrypoints returned help.
- All eight files rejected `python -O` before work. A separate check confirmed
  all eight also reject `python -OO`.
- `wedge5_spherical.json` reconstructed the frozen exact source byte-for-byte,
  SHA-256
  `28585450c46b4fbdc6f9eefec55b7e6d29f18fc527b654d00359ed425cb1f949`.
  All 13 off-origin shifts were zero.
- Existing output paths were rejected without changing their contents.
- The dry-run search command matched every final-stage setting, including
  `wedge5_step.json` as seed, ridge 0.03, two rounds, and bundled local
  spherical budget and far boxes. Zero rounds were rejected; no solver ran.
- A matrix requiring an off-origin shift was rejected before writing output.
- Both packages were copied to a temporary location and invoked from a
  different working directory. Canonical exports of the witness, exact PSD
  payload, and prerequisites matched the sibling proof package. Export-only
  correctly reported complete verification as false.
- Mutated spherical-budget bytes were rejected before creating output by both
  the exporter and the search wrapper. The uncertified numerical checkpoint
  was rejected as an input to the frozen complete-replay command.

Separate input checks confirmed byte identity of the warm start, final
checkpoint, exact source, fiber rows, far boxes, and approved spherical budget
with their originals. The frozen original relative provenance path was
retained. The dimension-five conservative spherical/Rankin/centroid cap table
is 10/9/8/5/5/5/3/2/1 with 132 bins; no newer sharp-angle budget is used.
All 17 files in the earlier centroid discovery package still matched its
frozen manifest, confirming that it was not modified.

## Complete copied proof replay

Ran `python export_and_check.py --output <new temporary directory>` against
the sibling `spherical_lens_verification_package`. The command generated all
three canonical payloads anew, copied the verifier scripts and frozen source
and budget, and ran the complete proof in the new output directory.

The replay ended with `FROZEN_EXPORT_AND_COMPLETE_VERIFICATION_PASS` and
`complete_verification_passed: true`. It used the hardened standalone far
verifier, which independently binds the payload, physical hash, dimension,
all nine cap thresholds, and geometry certificates before subdivision.
All copied script hashes were checked against the final frozen sibling proof
suite, whose `SHA256SUMS` hash is
`76c30a85ef155f8bb45691a53f133362c922209bb55da8359e8f581af65bcf29`.

- Frozen source, payload and boundary conventions; exact PSD and Schur
  factorization; radial/Fourier signs and zero Gaussian corrections: passed.
- Exact algebra, conservative transverse projection, embedded spherical LP
  signs, Rankin/centroid caps, tail budgets, candidate binding and explicit
  half-open count endpoints: passed.
- Pair verifier: 123 checked boxes, maximum depth 11, no failure.
- Cyclic verifier: 99 checked boxes, 44 passed, 6 excluded, maximum depth 12,
  no failure.
- Far verifier: 2,879 checked boxes, 1,011 passed, 430 excluded,
  113 relaxed boxes, maximum depth 22, no failure.

The final outward density enclosure was:

```text
[0.5099534974153451580358966431758354258465687576453685354517475202369317522156 +/- 7.46e-77]
```

The strict upper is **0.509953498** for unrestricted sphere packings, with
published near-neighbor count bounds retained as external mathematical inputs.

Exact SHA-256 identities of the copied proof scripts used in this replay:

```text
ab1990bedec83941063d76c87af2887519ab96d69c07afed3db2068cf3f81920  verify.py
8852da12f3930c7baa94c86d8ee1f748b999c1627e8cec7d97ae76117a42adf9  gaussian_guards.py
5e7aad6e82cbad749ccd41c706661937517126145befa2276a7150d503975e01  audit_algebra_exact.py
05b162191ea15021710e4410536311dffe01c43f122cd52b06980fae305e5b0f  construct_step_allowances.py
d51f5c51932e1646a4ac6da0745c38acd248bcd7f62ed22edc808f40dc159858  step_allowance_query.py
990f8149a66e29e8a78fc5d2b23232d8cc29690679a76992ca82186f960fc20c  verify_count_endpoints.py
c3489bbf1e1cf97ecb96deb659751773923948727ff87e7a0675899e0143c539  certify_count_pair_arb.py
fc62a7c1ba88b3310f5d9b292eaf3292429c848217363983994adcdfd6770f48  certify_cyclic_bernstein.py
7ae3e416f364d23601369aabb7b9adbfabc8b0cb3e017f634091bc3d85ba27f0  certify_general_lens_far.py
0933a774c0ec6fad329b05216f5a5e1ed220a63b6344fb6724ebd3fa827936e2  spherical_lens_caps.py
```

## Limits

The numerical optimizer and full exploratory history were not rerun. Only
the final bounded warm-start inputs and command are supplied. Fresh numerical
outputs can differ and remain uncertified until separately prepared,
completely verified, and reviewed. The successful complete replay above
applies only to the frozen exact source.
