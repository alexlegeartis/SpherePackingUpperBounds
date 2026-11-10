# Dimension 48

Current claim: density < 4.13385739314e-7 for unrestricted sphere packings. Read continuum_three_point/CYCLIC40_FULL48_PROOF.md and INDEPENDENT_CYCLIC40_FULL48_AUDIT.md. Other proof notes are historical dependencies, not additional registered claims.

Run from this directory with Python 3.11+ and python-flint 0.9.0 installed:

    python continuum_three_point/independent_cyclic40_angular.py
    python continuum_three_point/independent_cyclic40_domain.py
    python continuum_three_point/independent_cyclic40_radial.py

Do not use Python -O. No optimizer is required. The cited continuum theorem is an explicit mathematical input.
