# Stronger unrestricted 52D sphere-packing upper bound

**Every sphere packing in R^52 has density less than 9.6610185 × 10^-8.**

The degree-29 three-feature compact cyclic harmonic-matrix kernel uses radial powers 12, 18, 24, a degree-48 compact correction and degree-71 Gaussian/Fourier polynomials. No lattice or antipodal restriction is imposed. The full independent 512-bit geometry replay covers 100,505 nodes / 50,253 leaves. Exact PSD and origin conditions, independently reconstructed six Bernstein controls and scalar profile, angular edges, Fourier/tail signs, 629 annuli/caps and even-dimensional normalization pass. See the archived AUDIT.md and THEOREM_SCOPE.md for proof premises and independent/shared-code scope.

## Reproduce

Install `python-flint` and `sympy`, keeping assertions enabled:

```sh
python run_all.py
python run_all.py --full
```

The runner verifies the archive and every stored file before safe extraction. Full mode regenerates the radial annuli and independently replays the complete geometry trace. Quick mode checks exact bindings, PSD, angular controls/edges, Fourier/tail, budgets and normalization. The archive includes source, frozen data, trace, all seven completed audit logs and checksums. Standard positive-definite-kernel/Poisson theorems are mathematical inputs, not formally verified; shared scalar code was reviewed and rerun, not wholly independently reimplemented.

## Comparison

This improves the preserved repository certificate <9.7440293e-8 by **0.85191452%**. Relative to the [Afkhami-Jeddi et al. numerical spinless-bootstrap estimate](https://doi.org/10.1007/JHEP12%282020%29066) 1.0661233909836e-7 (N=466), the reduction is **9.3817978134%**. [Source dataset](https://hdl.handle.net/1721.1/125646). That external value is numerical, not a certified prior bound or exhaustive worldwide-priority guarantee.
