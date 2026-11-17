# Certified unrestricted sphere-packing upper bound in dimension 57

**Every sphere packing in R^57 has density less than 1.5718326 × 10^-8.**

**Comparison unresolved: this is not a SOTA claim.** The certificate establishes the stated inequality, but comparison with the strongest applicable prior bound has not been settled. No improvement percentage or record status is asserted.

The degree-29 two-feature compact cyclic harmonic-matrix kernel has exact PSD blocks and a certified origin extension. The independent 512-bit replay covers 121,409 nodes and 60,705 accepted leaves. Full angular signs, 12 edges, degree-71 Fourier/tail proofs, 621 annuli from contact 25 through cutoff 180, exact code-cap budgets and odd-dimensional Gamma/pi normalization pass. No lattice or antipodal restriction is assumed.

## Reproduce

Install `python-flint` and `sympy`; keep assertions enabled.

```sh
python run_all.py
python run_all.py --full
```

The runner verifies the archive digest, safely extracts it to a temporary directory, checks every file hash and executes the verifier. Quick mode checks exact PSD/bindings, angular corners, edges, Fourier/tail conditions, code caps, normalization and the polynomial-shift validation. Full mode also regenerates every radial annulus and independently replays the full three-dimensional partition.

[proof.tar.xz](proof.tar.xz) includes all exact inputs, proof source, partition trace, seven original independent audit logs and receipts. `STORED_FILES.json` lists content digests. Historical producer notes are clearly separated in provenance. Frozen proof inputs are byte-identical to the independent audit. The shared scalar algorithms were reviewed and rerun, not wholly independently reimplemented; standard kernel and Poisson theorems are mathematical inputs, not formally verified here.
