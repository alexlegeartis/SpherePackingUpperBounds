# Changes from the private code

To see every change: `diff -u original/verify_finite_count.py run.py`.
`multiplicity/finite_count_gain.py` is byte-identical to the private file.
No other private file was edited.

## `original/verify_finite_count.py` to `run.py`

Original path:
`work/d4/d4-11-resume-20261022/local-packing-energy/verify_finite_count.py`
(run 2026-10-22; its output is `data/census_floors.json`).

| Change | Reason |
|---|---|
| `PROOF` (the stage 2 output) is `--stage2 DIR` or `../stage2_origin/reference`, instead of the private `run-final/proof`. With `--stage2`, `DIR/RESULT.json` must say `PASS`. | Public chain of stages. |
| `bind()` takes a path relative to the private `source/` tree. The file is read from `../stage2_origin/source/`, or from `../../data/v7/` for the candidate data. The digest key is that path with the prefix used by the stage 2 output: `source/` for a fresh run, `work/d4/radial-stratified-origin-bound/source/` for `reference/`. | Same identity check against the parent acceptance; only the file locations differ. |
| The candidate value `new_beta` comes from a fresh run of `multiplicity/finite_count_gain.py` on the same stage 2 output, written to `DIR/finite_count_gain.json`, instead of the stored private `finite_count_gain.json` (640 KB). The script is unchanged; `run.py` sets its module globals `ROOT`, `SOURCE` and `OUT` after import. | Recompute instead of replay. |
| Output goes to `DIR/FINITE_COUNT_ACCEPTANCE.json` instead of the private package folder. Its `source` field is the input folder relative to the D4 folder. | Plumbing. |
| `main()` takes the output folder and candidate path and returns the report. | Plumbing. |
| New `compare_with_census()`: the fresh report must equal `data/census_floors.json` in every field except `source`, with the 124 floors, `new_beta`, `exact_raw_objective`, `exact_objective` and `strict_upper_bound` compared as exact rationals. | Ties the published data to this run. Adds checks; removes none. |
| New `cli()`: `--out`, `--stage2`, fail-closed `RESULT.json`, exit 0 only on `PASS`. | Public interface. |

Every mathematical check of the original is kept with the same code.

## `finite_count_gain.md` (documentation)

Original path: `work/d4/d4-11-resume-20261022/local-packing-energy/multiplicity/finite_count_gain.md`.
The written derivation of the finite-count inequality. Three edits, all
outside the mathematics: the opening note no longer says the result awaits
coordinator verification; the reference to the stage 2 driver points to
`../stage2_origin/original/verify.py` (same lines 171 to 186); the run command
with a private absolute path is replaced by the `run.py` command.
