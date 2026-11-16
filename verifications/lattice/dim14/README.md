# 14D lattice-packing upper certificate

Every lattice packing in R^14 has density **at most .034712755778**. This is not an unrestricted sphere-packing upper bound or an optimality claim.

From the repository root run `python verifications/lattice/replay.py --dimension 14`. This materializes hash-verified compressed inputs into a temporary directory, then executes `run_all.py`. Python 3.10+ is required; never run with `-O`.

The complete chain is:
1. Exactly verify the original Zubrilina prior certificate with the vendored integer Taylor verifier.
2. Prove the arithmetic trace-gap alternative: either rho ≤ 512, or every 14 distinct opposite pairs in one modulo-two class have total squared length at least 28+12/25. The geometric determinant and finite torsion computations are included with their analytic proofs.
3. Independently verify the new rational auxiliary polynomial, global Fourier positivity and spatial tail, every interval majorant, and every first-14 profile dual case. The low-density alternative is also checked below the reported bound.

See [the common method](../METHOD.md), the proof documents in this folder, and [the separate lattice comparison table](../../../LATTICE_RESULTS.md). The determinant step is attributed to Géraud-Stewart and Naccache, Appendix D, Lemma 32 of https://eprint.iacr.asia/2026/2826.pdf . Independent reconstructed checks are not a novelty claim for that lemma.

The theorem is an analytic argument with exact finite subcomputations, not a formal proof-assistant development. Numerical optimizer status is not used as proof.
