# Original-prior certificate verification

`taylor.py` and `bernstein.py` are unmodified upstream verification programs from Nina Zubrilina's `AstroNinja389/lattice-packing-bound` repository, pinned at commit `340a0e93bb95942c2b1b80092be43149504970eb`:

- https://github.com/AstroNinja389/lattice-packing-bound/blob/340a0e93bb95942c2b1b80092be43149504970eb/verify/taylor.py
- https://github.com/AstroNinja389/lattice-packing-bound/blob/340a0e93bb95942c2b1b80092be43149504970eb/verify/bernstein.py

The per-dimension `prior_certificate.json` files are original `certificates/n014.json`, `n015.json`, and `n018.json`. Source attribution: Nina Zubrilina, *Improved lattice-packing bounds in small dimensions*, https://arxiv.org/abs/2610.35830 . The original files and their claims are preserved; they are not new coefficients from this work.

Exact integer Taylor replays passed locally on 2026-11-16:

```
PASS [proven] n=14 K=61 density=0.0352392080 hermite=2.66846231
PASS [proven] n=15 K=61 density=0.0256118267 hermite=2.79036749
PASS [proven] n=18 K=61 density=0.0103704834 hermite=3.17829928
```

These prior upper bounds close the input assumptions used in the arithmetic trace-gap proofs. The new certificates and geometric arguments are verified separately. The 16D package includes its own earlier rational weighted-frame certificate and verifier, so it does not rely on an unchecked source-density metadata field.
