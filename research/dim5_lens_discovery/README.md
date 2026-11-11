# Five-dimensional lens discovery companion

Portable discovery and deterministic export companion for the independently
audited unrestricted sphere-packing bound **density < 0.509996147**. The proof,
theorem, external count-bound inputs, and complete verification live in the
sibling `lens_verification_package`. This directory is not a substitute for
that proof package. Numerical optimizer output is never automatically a bound.

## Setup and deterministic replay

Keep both package directories next to each other, or pass the proof directory
explicitly with `--verification-package`. Python 3.12.14 was used for the
recorded tests. The dependency versions are pinned in `requirements.txt`.

```sh
python -m pip install -r requirements.txt
python test_discovery.py
python reproduce_exact.py reconstructed
python export_and_check.py --candidate reconstructed/source_exact.json --output verified_replay
```

Run commands from this directory or invoke them by an absolute path from any
working directory. Choose new output paths. No remote action is performed.
Complete verification needs NumPy and python-flint, but no optimization solver.
The final numerical search additionally needs SciPy, CVXPY, and Clarabel.

`reproduce_exact.py` rebuilds `source_exact.json` from the saved numerical final
checkpoint `wedge5_step.json` and requires byte-for-byte equality with the frozen
source, including its original metadata and provenance. Its SHA-256 is:

```text
235920216afbd10b38be5efd3672e9d94b5ec5af8fc54aca3a63d77d9997721f
```

Preparation takes the exact rational average of each pair of finite-decimal
matrix entries. It requires all off-origin blocks to be positive definite,
allows **no off-origin diagonal shifts**, and sets only the isolated `A[0,0]`
to the exact Schur value plus `1/10^12`. All 13 recorded shifts are zero. No
kernel correction or nonzero reciprocal-Gaussian correction is introduced.
These are prerequisites, not the complete geometric verification.

`export_and_check.py` re-exports the frozen source, adjusts only the lens
metadata, and checks canonical JSON equality of `witness.json`,
`witness_psd_exact.json`, and `prerequisites.json` against the sibling proof
package. It copies the final proof scripts and local budget into the new output
directory and runs `verify.py` there. This checks the frozen source identity,
exact PSD and Schur identities, radial/Fourier conditions, centroid cap and tail
budgets, count thresholds, pair inequalities, cyclic inequalities, and full
lens-aware far domain. The accepted outward density enclosure is approximately
0.5099961464057082278, strictly below 0.509996147.

Use `--export-only` only to test canonical export equality. The resulting
report explicitly marks complete verification as false. The frozen verifier
deliberately rejects other candidates; a fresh numerical result needs its own
complete exact verification and review before any bound is asserted. Neither
preparation nor export alone confers that status. Python `-O`/`-OO` is rejected.

## Bounded final warm-start search

```sh
python run_search.py --output new_candidate.json --dry-run
python run_search.py --output new_candidate.json
```

This reruns only the **final bounded two-round warm-start stage**, not the full
exploratory history. The bundled `wedge5.json` is an uncertified intermediate
seed. The selected `wedge5_step.json` is the saved final numerical checkpoint.
Its solver status was `optimal_inaccurate`; exact downstream verification,
not the solver status or sampled residuals, establishes the frozen bound.
No fresh optimization was run as part of packaging. Solver outputs and timings
may vary even with these pinned versions and the fixed search RNG seed.

The wrapper supplies the following final-stage command, resolving every input
to its bundled local path:

```sh
python search_wedge_step.py \
  --dimension 5 --nodes 5 --angular 12 --degree 31 --ratio 2.25 \
  --subdivisions 2 --counts --schur --floor 1e-9 --ridge 0.2 \
  --margin 2 --prox-floor 0.0001 --boundary --exclusion 1.0 \
  --wedge-radius 1.08 --far-boxes far_boxes5_v2.json \
  --step-budget step_allowances_005.json --rounds 2 \
  --seed wedge5.json --output new_candidate.json
```

The original `step_budget` relative path inside the frozen numerical and exact
payloads is provenance, not a runtime lookup. It is retained byte-for-byte to
preserve the source hash. Runtime CLIs use the bundled `step_allowances_005.json`.
The proof witness uses `lens_budget: "step_allowances_005.json"`.

## Inputs and provenance

- `search_wedge_step.py`, `cyclic_fiber_rows_5.npz`, `wedge5.json`, and
  `wedge5_step.json` originate in `density_shell_features_20261111`.
- `source_exact.json` is the byte-identical `wedge5_step_proof/source_exact.json`.
- `step_allowances_005.json` is the original
  `long_edge_lens_cuts_20261111/step_allowances_005.json`, SHA-256
  `cd4af8aaac2b6d2b1f0559a843953ee0a9132543fb371a3fad4ab9d14e0a44e2`.
  It contains the old **centroid** caps 83/8/5/3/2/1, short radius 27/25, and
  104 dimension-five bins of width 1/200. No newer spherical budget is used.
- `far_boxes5_v2.json`, `prepare_candidate.py`, `export_low.py`, and
  `gaussian_guards.py` originate in `density_low_structural_20261111`.
- Portable-source changes: the search refuses optimized Python; preparation
  refuses off-origin repairs and output overwrites; the exporter has an
  argparse CLI. The search and exporter mathematics are unchanged.
- `test_discovery.py` runs deterministic compile/CLI, exact reproduction,
  portable export, search-command, no-overwrite, no-shift, candidate-rejection,
  budget-rejection, and optimized-Python rejection checks. It runs no solver
  and no full proof unless you separately invoke `export_and_check.py` without
  `--export-only`.

See `TEST_RESULTS.md` for the actual checks performed and `SHA256SUMS` for the
final package integrity manifest. Hashes provide identity, not mathematical
proof. The sibling proof package documents the theorem and dependencies.

For repository integration, the default suite resolves to `verifications/dim5/refinement_lens` when no sibling standalone package is present. The exporter selects `verify_primary.py`; the registry wrapper additionally runs the independent audit. The binary `cyclic_fiber_rows_5.npz` file is a frozen numerical search input, not a proof input.
