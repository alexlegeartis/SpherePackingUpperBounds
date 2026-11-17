# Unrestricted sphere-packing upper bound in dimension 54

**Every sphere packing in R^54 has density less than 4.8695034 × 10^-8.**

This degree-29 compact cyclic harmonic-matrix certificate uses two radial features, exact rational PSD blocks and a repaired degree-71 Fourier polynomial. It assumes neither lattice structure nor antipodality nor mod-2 restrictions.

The complete compact triangle domain was independently replayed at 512-bit precision: 95,299 nodes and 47,650 accepted leaves. All angular signs, the infinite scalar tail, 905 contiguous radial annuli from contact coordinate 24 through cutoff 250, exact code-cap budgets and safe density rounding were verified. Even-dimensional density normalization is exact rational arithmetic.

## Reproduce

Install `python-flint` and `sympy`; keep assertions enabled.

```sh
python run_all.py
python run_all.py --full
```

The runner validates the archive digest, safely extracts to a temporary directory, checks every unpacked file and executes the verifier. Quick mode checks the polynomial-composition validation, PSD/bindings, angular corners, all twelve edge intervals, Fourier/tail conditions, exact cap budgets and normalization. Full mode additionally regenerates all 905 annuli and independently replays the entire three-dimensional partition.

[proof.tar.xz](proof.tar.xz) contains the complete proof source, exact inputs, trace, seven independent audit logs, AUDIT.md and completion receipts. `STORED_FILES.json` lists their digests. All frozen audited inputs are byte-identical to the independent audit. The replayer uses independently reviewed outward-rounded Arb polynomial composition for shifts; a separate exact rational validation is included. The scalar arithmetic algorithms were reviewed and rerun, not wholly independently reimplemented.

## Comparison and scope

The [Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini (2020) numerical spinless-bootstrap benchmark](https://doi.org/10.1007/JHEP12%282020%29066), truncation N=484, is approximately 5.1522853741015e-8; this result lies **5.4884765414%** below that numerical value. [Original dataset](https://hdl.handle.net/1721.1/125646). This is a comparison with a numerical estimate, **not a certified prior sphere-packing upper bound**. It does not certify that optimizer or establish exhaustive worldwide priority.
