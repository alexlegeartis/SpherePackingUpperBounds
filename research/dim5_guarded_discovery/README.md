# 5D discovery and exact export companion

This portable warm-start package accompanies the independently audited strict
upper **0.510129239**. Numerical output is not a bound.

    python -m pip install -r requirements.txt
    python run_search.py --output new_candidate.json
    python reproduce_exact.py reconstructed
    python export_and_check.py --output replay --verification-package ../../verifications/dim5/refinement_guarded

The final numerical stage starts from the supplied complete checkpoint and
saved separation cuts. Earlier joint radial/kernel optimization is described
in the research report; the minimal package does not replay the whole
exploratory chain from empty initialization. The included generic search also
supports joint radial coefficients. It builds exact-derived cyclic boundary
rows when needed, preserves the two common local far boxes, and adds the
dimension-specific far boxes stored here.

The selected search uses degree 31, five radial features, angular degree 12,
R=3/2, a fixed origin Schur relation, PSD floor 1e-8, three-times-larger sign
margins and a regularized objective. Fresh optimization may vary across
solver or BLAS environments. It is not certificate evidence.

The unguarded numerical source is INVALID. Exact preparation is followed by the supplied annulus-constant PSD kernel guard of weight 1e-8. The guard keeps the pair correction identically zero. Its exact polynomial and construction are included; the final proof checks all conditions directly.

reproduce_exact.py deterministically rebuilds the frozen rational source from
the saved numerical checkpoint and checks SHA-256
79881d16bc8d837d13222c0d3c5ca6d47264e031e09cd46738275198d4da85b1. The final source, not any earlier checkpoint, must be used for proof.
No off-origin regularization repair shift was needed.

Export without the verification-package option proves prerequisites only.
The full command recomputes every domain and the strict upper in the separate
proof package. Scripts have no absolute-path or repository dependency after
installing requirements. Never use python -O for exact checks.

In the integrated repository the portable source driver is `verify_primary.py`; export uses that driver instead of the registry wrapper. The full registered wrapper additionally replays the independent audits.
