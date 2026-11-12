# Independent 10D audit

These scripts do not call the primary verifiers. They read frozen copies of
the same mathematical inputs, with aggregate linkage checks.

- `verify_geometry_independent.py` adapts the separately implemented published
  analytic-gradient audit to dimension 10. Gegenbauer coefficients are
  reconstructed by closed form instead of the primary recurrence. Exact
  binary-tree reconstruction checks coverage, rather than relying on area
  alone. Analytic interval gradients verify all 15,220 accepted rectangles.
- `verify_radial_independent.py` uses the direct monomial formula
  Fourier(X^j exp(-X)) = j! L_j^4(Y) exp(-Y). Exact rational Sturm sequences
  establish the two global sign rays, without complex-root isolation.
  A first-derivative mean-value enclosure covers the pair condition in
  45 intervals, rather than the primary second-derivative Taylor bound.
- Both primary and independent computations use rigorous python-flint/Arb
  arithmetic; this is implementation independence, not an independent
  arithmetic library or machine-checked formalization of the analytic theorem.

Run the parent `verify_all.py` for all primary, independent, linkage, and
mutation tests. The explicit finite-sum theorem is reviewed in `../PROOF.md`.
