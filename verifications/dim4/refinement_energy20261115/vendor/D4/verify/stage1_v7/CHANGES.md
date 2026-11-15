# Changes to original files

Every file under `record/` is a byte-for-byte copy of the file at the same
relative path in the private source tree `work/d4/six-graph-v7-bound/source/record/`.
No file under `record/` was edited. Many of them are pinned by SHA-256 inside
the verifiers, so an edit would also break those bindings.

The four V7 data files are not copied into `record/`. `run.py` copies them
from `../../data/v7/` into the working record at
`current/candidates/round9_six_graphs_degree18/` after checking their
SHA-256 against `MANIFEST.json`. They are byte-identical to the originals.

The changes below are path plumbing only. None changes a polynomial, a
bound, a sign test, an acceptance test or a status string.

## `checks/check_v7_inputs.py`

- Original: `work/d4/six-graph-v7-bound/source/checks/check_v7_inputs.py`.
- Change: the expected SHA-256 of `verify_scalar.py` was read from the
  original `record/MANIFEST.json`. It is now the constant
  `SCALAR_VERIFIER_SHA256 = 069596bf...2d61e7`, the same value that manifest
  lists and that `audit/compose_six_graph_full_origin_credit.py` pins as
  `SCALAR_VERIFIER_SHA256`.
- Reason: that manifest lists all 2,957 files of the research archive, most
  of which are not part of this proof and are not exported.

## `checks/audit_scalar_origin.py`

- Original: `work/d4/2026-10-19-d4-v7-evidence-intake/origin/audit_scalar_origin.py`.
- Change: the hardcoded absolute project path and bundle path are replaced by
  a required `--record DIR` argument. The origin-credit module it imports for
  two negative controls is loaded from `checks/check_origin_credit.py` next to
  it instead of `evidence_review/origin_audit/check_origin_credit.py`.
- Reason: the original only ran inside the private project tree.

## `checks/replay_graphs.py`

- Original: `work/d4/2026-10-19-d4-v7-evidence-intake/graph/replay_graphs.py`
  (the graph reviewer's driver, status `GRAPH_LOCAL_ACCEPTANCE_PASS`).
- Change: the project root, source bundle, output folder and interpreter
  were fixed paths (`/Users/...`, `parents[4]`, a hardcoded Python). They are
  now `--record DIR` and `--output-dir DIR`, the current interpreter, and
  subprocesses run with the record as working directory. The check that the
  process runs from the project root and that `evidence_review/` exists is
  removed. The independent moderate checker is loaded from `checks/` instead
  of `evidence_review/graph_audit/`. Two report fields name the record
  instead of the source bundle.
- Reason: the original only ran inside the private project tree.

## Unchanged copies

- `checks/check_origin_credit.py` from
  `work/d4/six-graph-v7-bound/source/evidence_review/origin_audit/check_origin_credit.py`.
  Only its exact helpers (`sign_certificate`, `exact`) are used, by the
  negative controls in `audit_scalar_origin.py`. Its own `audit()` entry point
  is not run, because it compares against the historical monotone
  composition record, which this stage replaces with a fresh composition.
- `checks/verify_moderate_independent.py` from
  `work/d4/six-graph-v7-bound/source/evidence_review/graph_audit/verify_moderate_independent.py`.

## New files

- `checks/check_broad19_rays.py`: shows that the 5 broad19 nonorigin PSD
  rays of `origin.json` are integer Gram rays of `kernel.json`.
- `COMPONENTS.md` and `proofs/`: the component table and byte copies of the
  proof notes it cites (same relative paths as in the private record).
- `record/geometry/audit_six_graphs_small_cap_bounds_helper.py` and
  `record/geometry/SIX_GRAPHS_SMALL_CAP_BOUNDS_HELPER_CONSTRUCTION.json`:
  unchanged byte copies, read by `checks/replay_graphs.py`. The construction
  file holds old `/workspace` paths; `replay_graphs.py` rebinds them to the
  working record in memory, as the original driver did.

- `run.py`: the driver. It does what `source/replay.py --mode full` does
  (copy the record, run the coordinator, run the full composer) and adds the
  independent checks and a fail-closed final gate.
- `checks/quick_parent_checks.py`: quick mode only. It runs the two parent
  checks that the coordinator skips with `--stage pair` (spherical count
  certificates and scalar Sturm signs), using the original modules.
- `MANIFEST.json`, `README.md`, `CHANGES.md`, `requirements.txt`.

## Compiler

The native helpers `current/fast_split.py` and `current/factor_correlated.py`
call the program `g++`. `run.py` does not edit them. It puts a two-line shell
script named `g++` first on `PATH` that runs `clang++` (or `--cxx`). On macOS
`/usr/bin/g++` is already Apple clang, so this changes nothing there; on
Linux it selects clang++ instead of GCC.

## Hash bindings to files that are not exported

None were removed. The coordinator and composer bind only files they open,
and every opened file is exported (see `MANIFEST.json` and the read-set
section of `README.md`).

## Publication recovery

Restored the missing pointed-cap written derivation as `proofs/current/base/new_work/audit/POINTED_VERIFICATION.md` and registered it in the manifest. The manifest directory scan now ignores Finder `.DS_Store` metadata only; all other unlisted files still fail. This fixes the observed clean-copy preflight failure caused by OS metadata, without changing any mathematical check or acceptance inequality.
