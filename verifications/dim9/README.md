# 9D signed-profile upper-density certificate

Strict bound: **0.190941549**. Primary verification and separate geometry,
radial and normalization audits pass. Run `python verify_all.py` with
Python3.11+ and python-flint0.9.0. Never use `-O` or `-OO`. No optimizer or
network is required for replay.

The key idea permits signed hinge coefficients while keeping v*phi'(v)
nonnegative and nondecreasing. This broadens the earlier profile family
without losing all-scale reduction. See [proof and ideas](PROOF.md).

The average-neighbor count373 is an explicit external published theorem;
its original SDP is not replayed here. See EXTERNAL_COUNT.md. The package
is a fresh reconstruction, stronger than the lost0.190968050 certificate,
not byte-identical recovery and not an optimality claim.

Exact export from refit9.json is provided by export_exact9.py; profile and
angular input coefficients are frozen in geometry/. The joint optimizer and frozen discovery checkpoint are in
[../../research/signed_profile9_reconstruction/](../../research/signed_profile9_reconstruction/).

`SHA256SUMS` records the release snapshot. Rerunning checks may rewrite timing-bearing output reports; coefficient identities are checked by the replay itself.
