> Reproduction from the repository root: `python verifications/lattice/replay.py --dimension 18`. Large JSON inputs are compressed; the replay checks decompressed hashes in a temporary copy before running the proof programs below.

#18D exact certificate replay

Certified lattice packing density upper bound: **0.010343803448**.

Run `python run_all.py` from any directory. Python3 is sufficient; the vendored Taylor prior verifier and the new independent verifiers use only the Python standard library. Do not use `python -O`.

The replay checks, in order:
1. The original18D source certificate, establishing the prior density bound used to deduce rho<33100.
2. Exact determinant thresholds, with the dimension-uniform analytic argument in `DETERMINANT_PROOF.md`.
3. All small-torsion coordinate histograms.
4. All remaining prime/order25 cases, rank-two torsion and index arithmetic, including an exact Machin bound for pi. The geometric/group argument is in `TRACE_THEOREM.md`.
5. The new certificate: independent Laguerre reconstruction,960 interval inequalities,382 complete profile cutoff cases, global Fourier/tail signs, exact objective and the rho<=32768 alternative.

`certificate.json` is an unchanged copy of the audited candidate. Its SHA256 is
`54449bd3a1ec58acc56218dfcb4c7c944462bd201ff8b86fd41a7e38c57fb0a7`.

Global polynomial signs are proved by rational Bernstein subdivision after compactifying the half-line. This does not trust numerical root sampling or optimizer feasibility. The theorem remains an analytic proof with finite exact subcomputations, rather than a formal proof-assistant development.
