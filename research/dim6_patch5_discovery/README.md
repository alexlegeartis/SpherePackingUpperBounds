# Direct follow-on search and export

This is the portable companion for the separately audited strict bound
**0.409128533**. It does not modify the approved 99/100 mixture package.

    python -m pip install -r requirements.txt
    python run_search.py --output new_candidate.json
    python reproduce_exact.py prepared
    python export_and_check.py --output replay --verification-package ../../verifications/dim6/refinement_patch5

run_search.py builds its exact-derived cyclic-boundary discovery rows when
needed, then runs the bounded warm start from the supplied numerical source.
It preserves complete antipodal rows, three complete cyclic boundary fibers
and sixteen added far boxes. Numerical output is not a certificate and may
vary across environments. The saved original run stopped after one round,
with no sampled violations, but complete exact verification was still required.

reproduce_exact.py instead uses the frozen numerical checkpoint, symmetrizes
its matrices exactly, and sets the isolated origin to its exact Schur value
plus 1/10^12. No off-origin diagonal shift was required. It reproduces the
frozen source byte-for-byte and checks its hash. Export/replay then recomputes
every complete proof condition. No Gaussian repair is needed.

The verifier's conclusions do not depend on trusting the numerical Bernstein
rows, the solver, or sampled constraints. Standard analytic facts and the
published near-neighbor count inputs remain external dependencies. No remote
action or publication is performed by these scripts.

The repository export wrapper selects `verify_primary.py` from the integrated proof package; `verify.py` additionally binds the registry and runs independent checks.
