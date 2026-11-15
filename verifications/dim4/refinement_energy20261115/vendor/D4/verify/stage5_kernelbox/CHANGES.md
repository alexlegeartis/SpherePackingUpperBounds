# Changes to the original code

The scripts in `src/` are copies of the project's working files. Only file locations changed. No formula, tolerance, subdivision rule, gate or pass condition changed, except the two input substitutions marked "input" below. `paths.py` is new. It holds every input and output location, relative to the package.

Original locations are in the private project tree. Abbreviations: `KB` = `work/d4/d4-ideas-fleet6-20261022/local-packing-energy/d4-12-kernelbox`, `LPE` = `work/d4/d4-ideas-fleet6-20261022/local-packing-energy`.

| Package file | Original path | Change | Reason |
|---|---|---|---|
| `src/kformula.py` | `KB/kformula.py` | Removed `assert (ROOT / "CLAUDE.md").exists()` and `ROOT = Path.cwd()`. Read the V7 data from `data/v7/` and `ACCEPTANCE.json` from the stage 2 folder through `paths.py`. Docstring run line. | Run from any folder, no project files |
| `src/kcert.py` | `KB/kcert.py` | `sys.path` entry is the script folder. Docstring run line. | Same |
| `src/kbox_run.py` | `KB/kbox_run.py` | `sys.path` entry is the script folder. The project board stop control (`board.py control`) now runs only if `--participant` is given; the default is empty. | No board in the public package |
| `src/check_kernelbox.py` | `KB/check_kernelbox.py` | Script folder from `__file__`. Default trees from `data/kernel-box/`, default output in `paths.OUT`; the output folder is created. Docstring run line. | Same |
| `src/validate_kformula.py` | `KB/validate_kformula.py` | `sys.path` entry is the script folder (it held `kformula.py` and `kernel_poly.py` in two folders before). Output to `paths.OUT/validate_kformula.json`. | Same |
| `src/ledger.py` | `KB/ledger.py` | Census floors from `paths.CENSUS` (default `data/census_floors.json`, byte-identical to the original `FINITE_COUNT_ACCEPTANCE.json`). | Same |
| `src/exact_values.py` | `KB/exact_values.py` | Census floors from `paths.CENSUS`; output to `paths.OUT`. Input: the assertion that U_raw equals the value in the project file `shared_budget.json` now compares with `exact_source_objective` of the stage 2 `ACCEPTANCE.json`. | `shared_budget.json` is not public; the stage 2 file holds the same U_raw |
| `src/hrlemma_check.py` | `LPE/d4-12-hrlemma/hrlemma_check.py` | Removed the `CLAUDE.md` assertion. Kernel breaks from `data/v7/kernel.json`; output to `paths.OUT`. Input: check (5) read the run receipts `kbox_j1.status.json` and `kbox_j5.status.json`; it now reads the same fields (`status`, `lam`, `a`, `j`, and `failures == []` in place of `nfail == 0` and `state == "done"`) from the `summary` block of `data/kernel-box/kbox_j*.json`. | Receipts are not public. Check (5) is a presence check only; gate 4 of `check_kernelbox.py` replays every leaf |
| `src/coord_check_kernelbox.py` | `KB/coord_check_kernelbox.py` | Loads `kernel_poly.py` from the script folder. Output to `paths.OUT/coord_check_kernelbox.json` instead of a fixed path in another session's scratch folder. | Same |
| `src/controls.py` | `KB/controls.py` | Script folder from `__file__`; tree from `data/kernel-box/`; controls folder in `paths.OUT`. | Same |
| `src/kernel_poly.py` | `LPE/d4-12-cyclicslack/kernel_poly.py` | Removed the `CLAUDE.md` assertion. Loads `raw_v7.py` from the script folder. | Same |
| `src/raw_v7.py` | `work/d4/d4-11-resume-20261022/local-packing-energy/accounting/raw_v7.py` | `PROJECT = parents[5]` and the three project folders replaced by the script folder (helpers) and `data/v7/` (data). `ACCEPTANCE.json` through `paths.py`. Docstring run line. | Same |
| `src/energy_d4_gram.py` | `methods/families/local-packing-energy/checks/d4-target625/energy_d4_gram.py` | `PROJECT` is the D4 package root and `SOURCE` is `data/v7`. Only `raw_v7.py` imports this file, and it uses helper functions that do not read these paths. | Same |
| `src/pointed_graph_energy.py` | same folder as above | `SOURCE` is `data/v7`. Same remark. | Same |
| `src/independent_antipodal_allocation_exact.py` | `work/d4/radial-stratified-origin-bound/source/current/base/new_work/pointed/` | none (byte-identical) | |
| `src/two_capacity_antipodal_allocation_exact.py` | same folder | none (byte-identical) | |

The data files are not copied. `data/v7/*.json` are byte-identical to the project's `round9_six_graphs_degree18/*.json`, `data/census_floors.json` to `FINITE_COUNT_ACCEPTANCE.json`, `data/kernel-box/*.json` to the project's `kbox_j1.json`, `kbox_j5.json` and `ledger.json`, and `verify/stage2_origin/reference/ACCEPTANCE.json` to the project's `run-final/proof/ACCEPTANCE.json` (checked with `cmp` on 2026-11-13).

New files: `paths.py`, `run.py`, `negative_controls.py`, `README.md`, this file, `MANIFEST.json`, `requirements.txt` and the three notes in `proofs/`. The notes are edited from the project's `kernelbox.md`, `hrlemma.md` and the assembly handoff, with internal board references removed.

## Publication recovery: negative-control isolation

`negative_controls.py` now requires a new work directory, never removes an existing destination, limits workers to 1-4, and excludes `.lake`/Git/OS caches when copying the package. This prevents copying multi-gigabyte Mathlib caches after a Lean build. The altered certificate tests and their acceptance conditions are unchanged.
