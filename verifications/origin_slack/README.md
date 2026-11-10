# Verified origin-coordinate refinements

This package applies exact Schur-complement slack removal to five existing continuum witnesses. It changes only K(0,0), preserving every pair, cyclic, far, annular, Fourier and tail condition. All earlier source witnesses remain unchanged.

Strict outward density upper bounds:

- 4D: **<0.633428734445**
- 5D: **<0.510501094478**
- 19D: **<0.009618343139**
- 20D: **<0.006979497272**
- 21D: **<0.005149549312**

## Reproduce

From the repository root, install requirements.txt and run:

    python verifications/origin_slack/verify.py

This replays all eight original continuum verification paths and then performs independent Fraction LDL, exact Schur solves and residual factorization, and 256-bit Arb objective reconstruction with strict outward-rounding checks. Do not use Python -O. The algebra-only stage is verify_algebra.py; running it alone does not replay the inherited continuum proofs.

## Audit scope

Source witnesses are pinned to commit 85e395353a937369e5cff6400f66337967812270. The 4D source is refinement_p6, 5D refinement_p4; 19D–21D use their original certificate.json. Exact source-input SHA256 hashes are recorded in source_inputs.json. The independent implementation does not use the producer's SymPy inversion. Its exact alpha and normalization are checked against schur_records.json; its objective is reconstructed from the original witnesses, not rounded previous bounds.

All eight source replays passed in the isolated audit; logs are retained. The external published neighbor-count bounds in dimensions 4 and 5 remain mathematical inputs, as documented in the original proof packages. This is ordinary PSD slack removal, not a new multipoint theorem or a claim of worldwide novelty. See PROOF.md for the full implication.

schur_records.json preserves the earlier conservative research outputs and their historical replay-status fields. The sharper bounds above are directly asserted by the independent algebra verifier after successful source replay. The result registry points here.
