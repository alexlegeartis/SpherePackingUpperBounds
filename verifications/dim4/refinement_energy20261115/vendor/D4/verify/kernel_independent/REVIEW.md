# Independent-kernel review, 2026-11-13

The final implementation was rerun after input-validation hardening. The
comparison and failure-control processes both exited **0**; atomic terminal
receipts recorded successful child outcomes. The owning process groups were
reaped, and a final process-argument check found no surviving owned job.

## Reproduction

From the D4 folder, with its declared NumPy/python-flint dependencies and a
C++20 compiler installed:

```sh
python3 verify/kernel_independent/run.py --out NEW_DIRECTORY --workers 1
python3 verify/kernel_independent/negative_controls.py --out ANOTHER_NEW_DIRECTORY
```

The actual run used `/Users/kevin/miniforge3/bin/python3` (Python 3.14.6),
NumPy 2.5.2, python-flint 0.9.0, on macOS 27.0.1 arm64. It ran under an
operational resource/termination guard, which is not needed to interpret or
reproduce the checker. The full comparison took **622.8 seconds**. Guarded
aggregate RSS peaks were 264,650,752 bytes for the full comparison and
300,875,776 bytes for the failure controls.

## Recorded evidence

* [RESULT.json](evidence/RESULT.json): PASS for 84 full and 28 exterior cases,
  1,059,688 nonzero rational coefficients, seven pair pieces, exact origin
  and raw-objective equality, 775 exact component-parameter equations and
  49 radial basis polynomial identities. The raw objective agrees exactly
  with stage 1 and the packaged census; its displayed approximation is
  0.630000291443570.
* [CONTROLS.json](evidence/CONTROLS.json): PASS. Changing one spherical weight
  in the independent implementation and changing a Gram-ray weight only in
  its input both caused exact coefficient/pair mismatches. Ten malformed
  data variants and four malformed captured tensors were rejected. The
  unchanged three-case subset had all compared objects equal but correctly
  failed the full-coverage gate.
* [TERMINATION.json](evidence/TERMINATION.json): both final child outcomes
  were exit 0, with no wrapper failure. Its paths/PIDs identify the original
  run, not additional runtime dependencies.

The 1,397 boundary checks, 3,728 correction sample points and 238 internal
consistency samples passed. Samples are diagnostics and are **not** the
continuum proof. The polynomial equalities are exact coefficient identities.
The nonpolynomial equality argument is
[COMPONENT_IDENTITIES.md](COMPONENT_IDENTITIES.md), using exact parameter
matching and explicitly treating closed/strict boundaries, including
irrational diagonal roots. The publication coordinator separately reviewed
its truncation inequality, vanishing arguments, allocation joins and gauge
cancellation against the displayed sources. This is source/written-proof
review plus exact rerun evidence, not a second independent continuum sign
certificate or an unconditional Lean proof of the concrete kernel.

## Dependency review

The default driver resolves all input paths from its location in D4. Its
inputs are its local implementation/manifest, `data/v7`'s four JSON files,
`data/census_floors.json`, and stage 1's manifest-listed record/data files.
A direct resolution audit checked **135 manifest-listed paths**, all inside
D4. It copies the stage 1 record into the new output directory and places
the four data files at their declared locations. Its extraction imports the
coordinator from that copy; the coordinator selects the included
`current/base`, so its historical-snapshot fallback is not used.

`kernel.py` imports no stage 1 or stage 5 code. The comparison side imports
stage 1 solely to obtain its actual assembled objects and validated component
parameters. The Python dependencies beyond the standard library are NumPy
and python-flint. The included stage 1 Python sources contain no imports of
scipy, sympy, cvxpy, mpmath, gmpy, sage, pandas, networkx, ortools or psutil.
A C++20 compiler builds the included backends in the output record. No
network fetch, parent-workspace proof, scratch Lean build or work-folder
job guard is a verification input. The optional `bonus_stage5.py` sample
script reads only the packaged stage 5 sources and is not part of this
proof comparison.

The original detailed outputs remain in the recovery work package. These
three compact evidence files preserve the final acceptance and failure
receipts inside the publication tree; the source/data suffice to regenerate
all detailed per-case outputs. This run was in the normal project environment;
the coordinator owns the separate clean-copy package validation.
