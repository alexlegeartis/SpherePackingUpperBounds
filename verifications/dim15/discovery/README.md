# Discovery provenance

Numerical discovery requires numpy, scipy, and python-flint. The full original angular checkpoint is ../a12_b16.json. Its exact binary coefficients, constant guard1/10000 and renormalization deterministically reproduce ../angular_exact.json. The included angular scaling reference is only a numerical coordinate choice.

The final radial coefficients are radial_cut12.json, rationalized with Fraction.limit_denominator(10**16); the origin multiplier receives exactly1000001/1000000. The proof verifier checks this export against the certificate. The final radial search used degree63, s15/8, a12,b16, physical cutoff12, normalized margin1e−8, and saved cuts. All root-guided rounds are preserved. A fresh numerical solve is not expected to reproduce identical floating coefficients.

Commands, from this directory: python screen.py --a 12 --b 16 --rounds 2; python search_radial.py --degree 63 --s 15/8 --a 12 --b 16 --tag radial_cut12 --margin 1e-8 --cut 12. Historical refinement scripts are retained as provenance; use the frozen parent verify.py for acceptance.
