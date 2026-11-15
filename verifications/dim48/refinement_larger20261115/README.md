# 48D larger-profile refinement, 15 November 2026

Strict density bound: **4.0843570355E-7**.
Preceding main: 4.13385739314e-7 at `6e7990802b8e045a78e3f54c99669d5fe9b5c437`.

The radial profile is `h(v)=v^(-a) exp(-b(v-1))`, with
`a=-192`, `b=216`, and minimum squared distance
`s=6`. The physical and Fourier polynomials have degree
111; the discovery correction space has degree 95.
The package uses no external neighbor counts.

See [the shared proof](../../profile_enrichment_20261114/PROOF.md).
Run `python verify.py` with python-flint installed and assertions enabled.
The wrapper checks the registry and source hashes, then executes all four
primary and independent angular/radial checks. Both paths share FLINT/Arb.
Discovery solver success alone does not certify a bound.
