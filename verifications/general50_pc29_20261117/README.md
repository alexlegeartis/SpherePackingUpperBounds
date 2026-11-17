# General sphere-packing upper bound in dimension 50

**Every sphere packing in R^50 has density less than 1.99548924 × 10^-7.**

This certificate uses a compact cyclic harmonic-matrix kernel of angular degree 29, with two radial features of degrees 12 and 24. It is unrestricted: no lattice, antipodal, or mod-2 assumption is used.

The kernel has exact rational PSD harmonic blocks and a certified origin block. The scalar Fourier transform is an explicitly nonnegative factored polynomial times a Gaussian. Full compact triangle-domain coverage, angular signs at the support boundary, all radial residual intervals, the infinite scalar tail, and exact spherical-code error budgets have been checked.

## Reproduction

Use Python with `python-flint` and `sympy` installed. Assertions must be enabled.

```sh
python verify.py
python verify.py --full
```

The first command checks exact input bindings, PSD and Fourier/tail checks, angular corner signs, code caps, density rounding, and the integrity of all files. The full command also independently replays every accepted triangle-domain box at 512-bit precision and regenerates all 633 radial annuli at 768-bit precision. Full replay takes several minutes.

`main_independent_result.json` records 32,905 partition nodes and 16,453 accepted leaves. All leaves were replayed at higher precision than the producer. `AUDIT.md` explains the theorem, normalization, interval arithmetic, and audit scope. The global interval machinery was reviewed and freshly executed, not entirely independently reimplemented.

Safe rational upper bound: `199548924/10^15`. The exact certificate's rational upper bound is strictly smaller. The previous public 50D certificate remains in its own verification directory; this package supersedes its numerical bound without deleting its provenance.
