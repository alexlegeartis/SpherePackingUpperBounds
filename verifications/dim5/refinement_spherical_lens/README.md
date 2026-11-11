# Five-dimensional spherical-lens density certificate

This frozen certificate proves unrestricted sphere-packing density **strictly
below 0.509953498**, subject to the published near-neighbor count inputs.
It preserves the older centroid package separately and unchanged.

## Replay

    python -m pip install -r requirements.txt
    sha256sum -c SHA256SUMS
    python verify_primary.py
    python independent/verify_all.py

Copy this whole directory anywhere. Inputs resolve relative to the scripts,
not the working directory. No solver or network is needed after dependency
installation. The complete verifier uses NumPy and python-flint. Budget and
spherical-code replay alone require only python-flint; SciPy is imported only
for optional construction, never verification. Python -O/-OO are rejected.
The directory must be writable for regenerated reports and cyclic arrays.

The fresh complete replay checks exact source and parameter binding, all 66
positive rational PSD pivots, direct Fourier identity, global Sturm signs,
Schur saturation, nonnegative count costs, strict outward density rounding,
all 132 five-dimensional budget inequalities, exact projection/code caps,
explicit half-open count endpoints and complete pair, compact cyclic and far
domains. The stored extra 6D/7D budgets are also replayed (396 total).
No stored success flag replaces a fresh mathematical check.

Fresh producer counts: pair 123, compact cyclic 99, far 2,879 boxes
(113 positive-allowance leaves and 430 valid geometry exclusions).
The density enclosure is centered at
0.5099534974153451580358966431758354258465687576453685354517475202369317522156.

The unchanged source_exact.json SHA-256 is
28585450c46b4fbdc6f9eefec55b7e6d29f18fc527b654d00359ed425cb1f949.
Only witness.json lens_budget metadata changes to the relative
step_allowances_spherical_005.json. The physical function is unchanged from
the direct-5D certificate; both Gaussian correction weights remain zero.

The early dimension-4 spherical caps are 10/9/8 at base lengths
1.50/1.55/1.60. This package uses the conservative transverse cosine bound
1-(1-D)/(2U), then transverse Rankin and centroid caps. It does not use the
newer sharper projection lemma. The allowance at exactly 1.50 is zero.
See PROOF.md and RANKIN_AND_SPHERICAL_PROOF.md for geometry and boundaries.

INDEPENDENT_AUDIT.md and independent_audit_results preserve the separate
completed audit. Pair, budgets, spherical polynomials and far were independently
reconstructed; the compact cyclic checker was source-reviewed and freshly
replayed rather than independently rewritten. The standalone audit is a
separately preserved companion. Historical theorem text is reference context,
not the specification of this nine-threshold input.

The original counts 45/46/47/48 remain external mathematical inputs from
https://arxiv.org/html/1206.2608 Table 2, with each radius reduced by 10^-6.
Their spherical SDP certificates are not rerun. The new lens code polynomials
are fully replayed here. No optimality, exhaustive literature priority, or
proof-assistant formalization is claimed.

SOURCE_MANIFEST.json records source and packaging provenance; SHA256SUMS
freezes all shipped inputs, code and documentation. Regenerated outputs with
runtime paths/timings/arrays are excluded from that immutable hash manifest.
TEST_RESULTS.json records fresh malformed-input, optimization, compilation,
and unrelated-working-directory full replay tests.

The repository wrapper binds the registered strict bound and frozen source/budget hashes, then runs the portable primary verifier (`verify_primary.py`) and the separate independent full audit. The earlier centroid lens proof is preserved in the sibling package.

Within the repository, `python verify.py` additionally checks RESULTS.json and the registered package path before running both complete chains. That registry wrapper requires the repository layout; the two commands above are portable standalone replays.
