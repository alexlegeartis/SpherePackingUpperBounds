# 33D combined-branch profile refinement, 15 November 2026

Strict density bound: **0.000089087196692**.
Preceding main: 0.00009009682341 at `187124532a7e1e33d157e59906b2b9eb83e053df`.

The radial profile is `h(v)=v^(-a) exp(-b(v-1))`, with
`a=-96`, `b=120`, and minimum squared distance
`s=33/8`. The physical and Fourier polynomials have degree
111; the discovery correction space has degree 95.
The package uses no external neighbor counts.

See [the shared proof](../../profile_enrichment_20261114/PROOF.md).
Run `python verify.py` with python-flint installed and assertions enabled.
The wrapper checks the registry and source hashes, then executes all four
primary and independent angular/radial checks. Both paths share FLINT/Arb.
Discovery solver success alone does not certify a bound.
