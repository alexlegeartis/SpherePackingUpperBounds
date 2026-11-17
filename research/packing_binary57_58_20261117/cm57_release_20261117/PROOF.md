# Interval-safe exhaustive CVP research verifier

No new packing minimum is certified merely by this tool's existence or by its small-instance tests. A node-capped run reports `NOT_CERTIFIED`. Only a completed `EXHAUSTIVE_NO_VECTOR` run can support a lower bound, and the implementation and data must receive independent review before a new packing claim.

## Exact input

The exporter constructs the rational LDL decomposition `G=L D L^T` and verifies every entry of that identity exactly. Each rational coefficient is enclosed between two binary64 numbers, with the enclosure checked by exact conversion of those numbers back to rational values. The `bits` input format serializes their uint64 bit patterns and the C++ reader copies those bits into doubles, avoiding decimal floating-point parsing assumptions.

For target `r/3`, the squared distance from an integer vector z is

`sum_i d_i (z_i − c_i)^2`, where `c_i=r_i/3 + sum_(j>i) L_ji (r_j/3−z_j)`.

The exporter also verifies an even integral Gram and `r G r^T ≡0 mod6`. Therefore every numerator `(r−3z)G(r−3z)^T` is divisible by6. In the MW80 experiment, excluding every numerator ≤84 proves distance squared at least90/9=10. The target coordinates are mapped exactly from the saved original coset family.

## Exhaustive branch bounds

The recursive routine processes coordinates from last to first. Its accumulated quantity is a lower bound on the squared contribution of already assigned coordinates. It computes an interval containing c_i and an upper bound on the remaining radius. Every integer in the resulting enclosing interval is visited. The new contribution uses a downward bound on distance to the center interval and a downward diagonal bound. Thus no vector with exact norm below the requested cutoff can be pruned.

All interval endpoints are widened with `nextafter`. The implementation requires IEEE binary64 nearest-rounding arithmetic, correctly rounded basic arithmetic and square root, and the documented compiler flags (`-ffp-contract=off -fno-fast-math -frounding-math`). No aggressive floating-point reassociation or fused contraction is permitted. Nonfinite arithmetic, unsupported coefficient/coordinate ranges, malformed input, or a node cap result in `NOT_CERTIFIED`, never in an absence claim.

At each reached leaf the original integer quadratic form is evaluated in signed128-bit arithmetic. Explicit bounds (`n≤110`, Gram entries at most2^30, residual coordinates at most2^22) keep the absolute sum below2^88, well within that accumulator. A reported witness is therefore checked exactly.

## Faster conservative center enclosure

The optional `fast` mode is not heuristic pruning. It requires target numerators in{0,1,2}, all lower-triangular coefficient intervals inside[−16,16], interval widths at most2^−47, and enumerated coordinates bounded by2^20. It computes the center with ordinary scalar operations and encloses it with radius2^−10.

Here is a deliberately loose error proof. Each exact coordinate difference has magnitude below2^21. A coefficient approximation contributes error at most2^−26 per term. Coordinate subtraction and multiplication rounding contribute at most2^−28 each per term. There are fewer than128 terms. Every intermediate accumulator has magnitude below2^32, so each addition contributes at most2^−21 and their total is below2^−14. Including target division rounding, the total center error is strictly below2^−13. The chosen radius2^−10 is therefore conservative by a factor exceeding8. The same subsequent interval recursion remains exhaustive.

Exact size reduction of the MW80 dual-reversed basis bounds the actual Gram–Schmidt coefficients by1/2 while preserving the exact GSO diagonal sequence. This improves arithmetic conditioning without changing the lattice or coset.

## Tests and current status

`test.py`, `test_fast.py`, and `test_bits.py` compare the ordinary interval, fast-center, and bit-input variants with exact brute force on40 small positive-definite integer instances each. These tests support implementation checking but are not substitutes for the argument above or independent review. The current80-dimensional runs are bounded research runs; inspect their explicit status before drawing any conclusion.
