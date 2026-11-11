# Discovery provenance

The original angular snapshot is ../source.json. It is reproduced exactly by rationalizing its stored binary coefficients, subtracting1/10000 from c0 and renormalizing. The stored input scaling reference affects discovery coordinates only.

Optional discovery command: python screen.py --dimension 12 --a 10 --b 16 --rounds 2. It generated the saved angular snapshot; optimizer replay may differ numerically.

The final radial model used normalized Laguerre degree63, s=12/8, a=10, b=16, cutoff12, normalized row margin1e-8, and the saved cut list. Command: python search_radial.py --dimension 12 --degree 63 --s 12/8 --a 10 --b 16 --cut 12 --margin 1e-8 --tag radial_cut12. Run from this directory.

Exact export uses Fraction.limit_denominator(10**16) on stored radial coefficients and origin cost, with exactly1000001/1000000 added to M. The parent verifier checks that export against the certificate, independently reconstructs the monomial Fourier transform, and proves every global gate. Histories record rejected intermediate candidates, not bounds.
