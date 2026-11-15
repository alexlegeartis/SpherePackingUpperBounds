# Release verification — 15 November 2026

The base proof was freshly replayed in full mode with Python 3.12.14,
NumPy 2.5.2, python-flint 0.9.0 and GCC, using two cyclic/box workers.
The run took **2266.2 seconds** and passed all seven stages. Its 277 source
files were checked byte for byte against the frozen vendored copy; see
`source_identity.json`.

The full replay covered all 112 cyclic branches, 165,503 count vectors,
the finite-count and graph refinements, the original box trees, both kernel
comparisons, **1,059,688 polynomial coefficients**, 25 ordered radial cells
and 25 endpoint cells. The unchanged upstream formalization boundary still
applies; this is a non-Lean computer-assisted proof.

The final extension invocation consumes that run's fresh census and stage-2
output. Its `base_replay_reused` flag is true because the base and extension
were run in separate invocations. The base was genuinely executed during
this release; it was not substituted by the upstream committed success
reports. `base_SUMMARY.json`, the seven stage reports and `base_driver.log`
record that run. `composition_RESULT.json` records the extension checks.

The new proof adds 308 interval leaves, ten exact central-sign leaves,
245 exact Gram range equations and the recomputed 124-group ledger.
Six malformed inputs were rejected; their log is retained. Exact source
reconstruction/export and the generated registry/comparison checks passed.
The historical origin algebra still passes with the new 4D registration.
Other dimensions' proofs were unchanged and were not rerun in this release.

These reports are historical evidence. Reproduce the theorem with the
default full command in [the package README](../README.md); reports alone
are not a replacement for running the verifier. Absolute filenames in the
reports are provenance strings and are not required certificate inputs.
