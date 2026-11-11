# Discovery and exact provenance for the improved 6D certificate

This portable companion documents the search leading to the verified strict
upper **0.409135399**. The proof package remains the authority for validity.
No numerical checkpoint, sampled objective or solver status is a bound.

## Deterministically rebuild the final rational source

    python -m pip install -r requirements.txt
    python reproduce_exact.py reconstructed

This uses the frozen selected numerical checkpoint, exact symmetrization and
Schur-origin preparation, and a rational 99/100 mixture with the previous
source. It reproduces the final source byte-for-byte and checks SHA-256
`41c60fae103355323abdfab8740c64b1089732bc0e0903fd27da845ed4259810`.
No optimizer is needed for that step. Every off-origin block was already
positive definite; the preparation made no diagonal shifts.

To export and replay the complete proof using the separate verifier package:

    python export_and_check.py --output replay --verification-package ../../verifications/dim6/refinement_mixture99

All outputs go into a new directory. This does not overwrite the proof
package. Without the verification-package argument, export only checks
prerequisites and is not a full certificate.

## Re-run numerical discovery

    python run_search.py joint --output new_joint.json
    python run_search.py refine --output new_refinement.json

The first command runs the three-round joint radial/kernel stage from the old
published source. The second restarts the four-round fixed-radial refinement
from its supplied frozen joint checkpoint. To chain a fresh first-stage
output, add `--seed new_joint.json` to the second command. Every new numerical
result requires new exact preparation and full verification.

The discovery script fixes the origin Schur relation, reduces artificial PSD
regularization to 1e-8, retains saved pair/far/cyclic cuts, and includes the
antipodal and two local far-region Bernstein constraints. The joint stage uses
degree-31 radial coefficients as variables. The refinement freezes them and
uses a regularized objective with three-times-larger sign margins. These
improve the search in the same analytic family, not a new density theorem.

The refinement's fourth checkpoint was selected because its complete cyclic
domain passed and its remaining far violations were small. It is still an
invalid unmixed candidate. The final rational mixture passes all exact
conditions directly and needs no Gaussian repair. Later unsuccessful mixture
trials and unrelated experiments are intentionally excluded from this minimal
companion.

The source's active mathematical paths match the successful search settings.
Numerical reruns are not guaranteed bitwise-identical across machines or
solver/BLAS versions. The frozen exact preparation is deterministic and tested
separately. Pinned package versions, settings and file hashes are recorded in
requirements.txt, provenance.json and SHA256SUMS. The scripts have no absolute
path or repository dependencies. Never use python -O for exact checks.

In this repository the standalone producer verifier is preserved as `verify_primary.py`; the exporter selects it instead of the registry wrapper. The registered wrapper also runs the independent audit.
