# 45D profile refinement, 14 November 2026

Strict density bound: **0.0000012656297919**.
Preceding main: 0.000001287489137 at `8a080d6a00909391a91ecbac097e4004c404c40d`.

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
