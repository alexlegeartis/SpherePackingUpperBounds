# Changes to the original code

None. The four files in `src/` are byte-identical copies of the project files (checked with `cmp` on 2026-11-13):

| Package file | Original path |
|---|---|
| `src/gap_check.py` | `work/d4/d4-allpackings-0625/broad-antipodal/gap_check.py` |
| `src/verify_annular_count_151.py` | `work/d4/graph-checkpoint-2026-10-19/source/geometry/verify_annular_count_151.py` |
| `src/verify_exterior_graphs.py` | `work/d4/graph-checkpoint-2026-10-19/source/geometry/verify_exterior_graphs.py` |
| `src/annular_count_151_certificate.json` | `work/d4/graph-checkpoint-2026-10-19/source/geometry/annular_count_151_certificate.json` |

`run.py`, `README.md`, `MANIFEST.json`, `requirements.txt`, this file and `proofs/graph_slack.md` are new. The proof note combines the project's `slack-1208.md`, `nine-point-gap.md`, the 151-count section of `ANTIPODAL_EXTERIOR_GRAPH_THEOREM.md` and the projection step of the 1358-edge graph proof, with internal board references removed. The release now also ships `src/verify_projected_eight.py`, which replays the projected Legendre certificate from `data/projected_eight.json` exactly. Its tests reject altered cutoffs, negative coefficients and JSON floats. The old cited nine-point theorem is retained for attribution only and is no longer used by the stage driver.

The project's `graph_recomposition.py` is not copied. It read a project file (`BASE_INPUT_TRANSFER.json`) that is not public. `run.py` makes the same check from `data/v7/origin.json` and `data/v7/graph_energies.json`: the graph charge 4769/25000000 and U_raw = objective + charge/8.
