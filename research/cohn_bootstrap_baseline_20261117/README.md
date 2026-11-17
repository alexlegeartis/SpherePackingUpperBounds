# Numerical spinless-bootstrap benchmarks, dimensions49–96

Source: N. Afkhami-Jeddi, H. Cohn, T. Hartman, D. de Laat and A. Tajdini, [High-dimensional sphere packing and the modular bootstrap, JHEP12(2020)066](https://doi.org/10.1007/JHEP12%282020%29066), [paper](https://arxiv.org/abs/2006.02560v2), [authors’ original MIT dataset](https://hdl.handle.net/1721.1/125646).

These are **published numerical spinless modular-bootstrap benchmarks**, converted to sphere-packing-density units. They are **not certified packing upper bounds or achieved packing densities**. The authors’ data contain optimized scaling dimensions and degeneracies, not globally sign-certified Cohn–Elkies auxiliary functions. The general same-radius lift from the minus-one Fourier eigenfunction problem to a Cohn–Elkies function is conjectural (paper Section2.3); positivity certification is a separate requirement. Do not label these entries rigorous inequalities or use convergence alone as an error certificate.

Each row selects the largest supplied truncation order N for integer c_total=n. The method solves the finite spinless-bootstrap equations numerically; it does not enumerate packings. Dataset c_total equals the Euclidean dimension, not half that dimension. The first gap Delta1 gives r²=2Delta1, and the conversion is (pi*Delta1/2)^(n/2)/Gamma(n/2+1). Displayed decimal digits are numerical approximations, not rigorous outward rounding. Successive-N changes in JSON/CSV are convergence diagnostics only.

[Complete49–96 table](../../NUMERICAL_BENCHMARKS49_96.md). Machine-readable values are in baseline49_96.json and baseline49_96.csv. This package contains no new proof claims or purported improvements over these benchmarks.

## Reproduction and provenance

The complete original dataminusone.txt was checked against the MIT-published MD5 and the SHA256 in SOURCE_MANIFEST.json. It is not duplicated in the repository. gap_history1_96.json records the first gap at every supplied order, allowing highest-N selection and8D/24D/48D controls. The current publication was independently re-extracted from the checksum-matched full dataset before table generation.

Run `python -B verify.py` with mpmath installed to recompute highest-N selection, density conversion at90 and120digits, controls, JSON/CSV agreement and the displayed table. To reproduce extraction from scratch, obtain the original file at the recorded source URL, put it beside extract.py and run `python extract.py`, then the verifier. No solver rerun or global positivity certificate is claimed.
