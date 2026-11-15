# Non-minimum-basis section screen — 16 November 2026

## Outcome

No packing improvement or candidate P48p section was found in this run.
A reproducible randomized discovery screen generated 50,000 abstract rank-five
even integral Gram forms with diagonal entries in {6,8,10,12,14,16,18} and
off-diagonal entries in [-7,7]. The fast screen retained forms whose minimum
over coefficient vectors in [-2,2]^5 was at least 6.

Only three forms with determinant below 1,800 survived that screen. Each was
then checked exactly: a finite coefficient box derived from
|v_i|^2 <= 6 (K^{-1})_ii contains every possible vector of norm below 6.
All three failed, with exact minima 4, 2 and 2, respectively.

This is a deliberately limited randomized screen, not an exhaustive search
over Gram forms. It does not establish a nonexistence theorem.

## Why this is not yet a packing result

Even a valid abstract Gram matrix would need an explicit primitive embedding
into the 48D ambient lattice, followed by exact verification of the section,
its dual translate cluster, and the resulting density. No such embedding was
found or tested here. The experiment therefore makes no SOTA claim.

## Reproduction

From the repository root:

```bash
python -m pip install numpy sympy
python research/nonminimal_sections_20261116/search.py
```

The seed and sample count are fixed in the script. Floating eigenvalues and
the bounded [-2,2]^5 screen are discovery filters only; the exact final
minimum check uses SymPy integer/rational arithmetic.
