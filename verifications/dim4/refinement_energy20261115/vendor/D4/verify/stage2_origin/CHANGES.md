# Changes from the private code

No file under `source/` was edited. All 32 files there are byte-identical to
the private tree, and `run.py` checks each one it binds against the
producer's `source/MANIFEST.json`. The only edited program is the driver.

To see every change: `diff -u original/verify.py run.py`.

## `original/verify.py` to `run.py`

Original path: `work/d4/radial-stratified-origin-bound/verify.py` (accepted run
of 2026-10-20, output in `reference/`).

| Change | Reason |
|---|---|
| Path constants. `PROJECT`, `ROOT`, `SOURCE`, `OC` point into a fresh working copy `DIR/tree/source`, built by the new `build_tree()` from `source/` and `../../data/v7/*.json`. | The original read a 393 MB private tree at fixed project paths. The original modules compute their own root from their file location, so they must run from one tree that holds both the code and the V7 data. |
| The stage 1 gate `accept_full_replay.check_run` on the private replay folder and receipt is replaced by `stage1_inherited()`. With `--stage1 DIR` it requires `RESULT.json` status `PASS`, mode `full`, 112 branches rerun, the same four V7 files in `DIR/record`, a completed full composition, and equal raw objectives. It then takes the old credited objective and threshold from `DIR/record/replay/v7_full_origin_composition.json`. Without `--stage1` it takes them from `reference/INHERITED_FULL_V7_ACCEPTANCE.json`. | Stage 1 is a separate public stage. The private gate checked a receipt with absolute paths of one machine. The values consumed (`objective`, `strict_upper_bound`, 112 branches) are the same. |
| The binding loop over `accepted['proof_file_sha256']` binds only the 11 entries this stage reads (`STAGE2_PROOF_FILES`): the four V7 data files, the three spherical count certificates and four count checker modules. Each digest still comes from the accepted record. | The other 109 entries are inputs of the stage 1 sign proof (interval engines, native libraries, graph certificates, older candidates). This stage never reads them; stage 1 checks them. |
| The comparison of the four V7 files with the private replay record moved into `stage1_inherited()` and runs only with `--stage1`. | Same check against the public stage 1 output. Without stage 1 output there is nothing to compare; the manifest still pins the bytes. |
| `bind(Path(__file__))` removed. | The driver is outside the working tree; `bind` rejects paths outside it. `EXPORT_MANIFEST.json` records the driver digest. |
| `counts.BASE = ROOT / HIST` (same subtree, constant renamed). | Unchanged path, shorter expression. |
| The inherited acceptance is saved as `STAGE1_INHERITED.json` instead of `INHERITED_FULL_V7_ACCEPTANCE.json`. | The content now records which stage 1 source was used. |
| `check()` no longer creates the output folder; `main()` creates it before building the tree. `check()` also returns the inherited values. | Plumbing. |
| New `compare_with_reference()`: every fresh report must equal its copy in `reference/` after removing `seconds` fields and `source_bindings`; every fresh source binding must have the same digest as the binding of the same file in the accepted private run. | Shows the reference read by later stages is what this code produces from these inputs. Adds checks; removes none. |
| New `main()`: `--out`, `--stage1`, fail-closed `RESULT.json`, exit 0 only on `PASS`. | Public interface. |

Every mathematical check of the original is kept with the same code and the
same order: projection, count geometry, radial planes, the spherical identity
of the two energies, the three angular proofs, the enumeration, the credit,
the objective and the outward rounding.

## Files copied without change

`source/` (32 files), `original/verify.py`, and `reference/` (placed by the
coordinator as a byte copy of the accepted output).

## `PROOF.md` (documentation)

Taken from the section "Mathematical composition" of the private review
`work/d4/radial-stratified-origin-bound/REVIEW.md` (2026-10-20). Links to
private files are replaced by references to the paper; section headings were
added; a short paragraph on the 124 groups was added. The mathematics is
unchanged. The review's sections on supplementary checks, retention and
history are not included.
