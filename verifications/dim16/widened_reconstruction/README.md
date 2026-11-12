# Reconstructed widened-shell 16D upper bound

Strict upper bound: **0.024900975**.

Run `python verify_all.py` with python-flint and NumPy. It replays independent geometry, exact radial algebra/rays/full pair coverage, exact spherical PSD/pair checks, independent polynomial reconstruction and integer Bernstein full-domain covering. All inputs are local; no optimizer or network needed. The exact count objective7724.885513367348... yields integer7724 at cap501/1000. The exact cap-to-shell and normalization gates are mandatory.

Run `python negative_controls.py` to test rejection of count input corruption, radial coefficient corruption and an invalid smaller integer count7723.

See PROOF.md. Newly reconstructed, not byte-identical to the lost original. Discovery provenance: degree63 radial LP over published log-hinge geometry, exact normalization and epsilon1e-6 repair; spherical per-factor-row nonnegative PSD reweighting of the retained Schoenberg/local Gram construction. The original finite-sum count theorem is unchanged. The frozen verification does not require any discovery script.

`INPUT_HASHES.json` pins the mathematical inputs. `SHA256SUMS` is a snapshot of the release files; rerunning verification rewrites timing-bearing audit reports, so those output hashes may change without changing the inputs.
