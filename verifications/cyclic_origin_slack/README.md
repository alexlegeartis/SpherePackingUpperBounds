# Saturated cyclic bounds in dimensions 18-20

Strict unrestricted density upper bounds:

- 18D: <0.013224473
- 19D: <0.009521000
- 20D: <0.006936042

    python verifications/cyclic_origin_slack/verify.py

Install root requirements.txt first and never use Python -O. This command replays all three complete unsaturated cyclic source packages, then exact Fraction Schur elimination/residual factorization and 384-bit Arb objectives. Algebra-only replay is verify_algebra.py and does not replace source domain checks.

All earlier positive-definite source certificates remain untouched. Saturation is a separate derived PSD kernel, and does not modify off-origin inequalities. See PROOF.md and source_inputs.json. The independent read-only implementation was first tested against producer witnesses, then rerun here against the semantically identical registered sources. No worldwide-priority claim is made.

A second independent standard-library Fraction implementation, verify_rational.py, reconstructs every saturated residual entry and proves the outward thresholds using rational Machin bounds for pi (squaring odd-dimensional densities). It is invoked after the source hashes and registry checks; it has no Arb dependency. Its source continuum dependencies are still mandatory.
