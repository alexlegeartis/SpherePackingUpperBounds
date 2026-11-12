# Independent 16D transfer audit

`geometry/verify.py` reconstructs normalized Gegenbauer coefficients by
closed form, checks exact Sturm signs and all-scale hypotheses, reconstructs
the binary partition tree and checks every rectangle with analytic interval
gradients. It does not call the primary geometry implementation.

`verify_transfer.py` proves C_new/B_new = C_old/B_old exactly and recomputes
the density with rational Machin enclosures for pi. It uses the separately
replayed published radial and pointwise-count certificates. These inputs are
not replaced by saved PASS metadata: the parent `verify_all.py` explicitly
runs the complete original primary and independent chains before this audit.

The conditional stronger radial/count work is not included here.
