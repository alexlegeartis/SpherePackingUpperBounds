# Five-dimensional spherical-lens discovery companion

Portable discovery and deterministic export companion for the independently
audited unrestricted sphere-packing bound **density < 0.509953498**. The
theorem, external count-bound inputs, exact spherical-cap certificates, and
complete verification live in the sibling `spherical_lens_verification_package`.
Numerical optimizer output is never automatically a bound. This package is
separate from the earlier centroid-lens package and does not modify it.

## Setup and deterministic replay

Keep this package next to `spherical_lens_verification_package`, or pass that
directory explicitly with `--verification-package`. Python 3.12.14 was used for
the recorded tests. Dependency versions are pinned in `requirements.txt`.

```sh
python -m pip install -r requirements.txt
python test_discovery.py
python reproduce_exact.py reconstructed
python export_and_check.py --candidate reconstructed/source_exact.json --output verified_replay
```

Commands work from this directory or by absolute script path from any working
directory. Choose new output paths. No remote action is performed. Complete
verification needs NumPy and python-flint, but no optimization solver. The
numerical search additionally needs SciPy, CVXPY, and Clarabel.

`reproduce_exact.py` rebuilds `source_exact.json` from the saved final numerical
checkpoint `wedge5_spherical.json` and requires byte-for-byte equality with the
frozen source, including its original provenance. Its SHA-256 is:

```text
28585450c46b4fbdc6f9eefec55b7e6d29f18fc527b654d00359ed425cb1f949
```

Preparation takes the exact rational average of finite-decimal matrix entries,
requires positive-definite off-origin blocks, allows **no off-origin diagonal
shifts**, and sets only the isolated `A[0,0]` to the exact Schur value plus
`1/10^12`. All 13 recorded shifts are zero. No kernel correction or nonzero
reciprocal-Gaussian correction is introduced. This is prerequisite preparation,
not complete geometric verification.

`export_and_check.py` re-exports the frozen source, changes only the lens
metadata, and checks canonical JSON equality of `witness.json`,
`witness_psd_exact.json`, and `prerequisites.json` against the sibling proof.
It copies the final proof scripts and bundled budget into the new output
directory and runs `verify.py` there. The proof checks source identity, exact
PSD and Schur identities, radial/Fourier signs, spherical LP and Rankin caps,
centroid caps, tail budgets, count thresholds, pair/cyclic inequalities, and
the entire lens-aware far domain. The outward density enclosure is
approximately 0.5099534974153451580, strictly below 0.509953498.

`--export-only` checks canonical export equality without complete verification;
its report explicitly sets complete verification to false. The frozen verifier
rejects any other candidate. A fresh numerical output needs its own complete
exact verification and review before a bound is asserted. Python `-O`/`-OO`
is rejected.

## Bounded final warm-start search

```sh
python run_search.py --output new_candidate.json --dry-run
python run_search.py --output new_candidate.json
```

This runs only the **final bounded two-round warm-start stage**, not the full
exploratory history. `wedge5_step.json` is the saved numerical seed from the
earlier centroid stage; `wedge5_spherical.json` is the selected final numerical
checkpoint. Numerical files carry no certification by themselves. The final
solver status was `optimal_inaccurate`; exact downstream verification, not
the solver status or sampled residuals, establishes the frozen bound.
No fresh optimizer was run during packaging. Results and timings may vary
even with pinned versions and the fixed search RNG seed.

The wrapper supplies the following final-stage command, using bundled local
input paths and checking the spherical budget hash before running:

```sh
python search_wedge_step.py \
  --dimension 5 --nodes 5 --angular 12 --degree 31 --ratio 2.25 \
  --subdivisions 2 --counts --schur --floor 1e-9 --ridge 0.03 \
  --margin 2 --prox-floor 0.0001 --boundary --exclusion 1.0 \
  --wedge-radius 1.08 --far-boxes far_boxes5_v2.json \
  --step-budget step_allowances_spherical_005.json --rounds 2 \
  --seed wedge5_step.json --output new_candidate.json
```

The original `step_budget` path and `wedge_caps` metadata inside the frozen
numerical and exact payloads are retained byte-for-byte as historical
provenance. Runtime allowances come from the bundled spherical budget, not
those historical cap labels. The proof witness uses the portable relative
path `lens_budget: "step_allowances_spherical_005.json"`.

## Spherical budget and provenance

The bundled budget is exactly
`long_edge_lens_cuts_20261111/step_allowances_spherical_005.json`, SHA-256:

```text
605a15c07f305a590bb1e0e2c80a9a88f3991023129bcbb5a88c103aa4552046
```

It uses the approved **conservative transverse projection, exact spherical LP,
Rankin, and centroid** geometry. It does not use the newer sharp angular
formula. The short radius is 27/25. Dimension five has 132 bins of width 1/200.
The nine cap thresholds are 3/2, 31/20, 8/5, 41/25, 17/10, 7/4, 9/5,
46/25, and 48/25, with caps 10/9/8/5/5/5/3/2/1 respectively.

- `search_wedge_step.py`, `cyclic_fiber_rows_5.npz`, `wedge5_step.json`, and
  `wedge5_spherical.json` originate in `density_shell_features_20261111`.
- `source_exact.json` is the byte-identical
  `wedge5_spherical_proof/source_exact.json`.
- `far_boxes5_v2.json`, `prepare_candidate.py`, `export_low.py`, and
  `gaussian_guards.py` originate in `density_low_structural_20261111`.
- Portable-source changes: search rejects optimized Python; preparation
  rejects off-origin shifts and output overwrites; exporter has an argparse
  CLI. Search and export mathematics are unchanged.
- `test_discovery.py` runs deterministic compilation/CLI, exact reproduction,
  portable canonical export, search-command, no-overwrite, no-shift,
  candidate-rejection, budget-rejection, and optimized-Python rejection checks.
  It runs no solver and does not run the complete proof.

See `TEST_RESULTS.md` for actual checks and `SHA256SUMS` for package integrity.
Hashes establish identity, not mathematical proof. Refer to the sibling proof
package for the mathematical argument, independent audit, and external inputs.

The integrated export path selects `verify_primary.py` rather than the registered wrapper, which additionally runs independent checks. A sibling standalone proof package remains supported. `cyclic_fiber_rows_5.npz` is a frozen numerical search input, not a proof input.
