# 13D audit verdict: PASS

The mathematical source is `../dim13_curved_ready_20261111/`. Only copies in this audit directory were written. The original certificate and angular data still match the frozen input hashes. The resulting bound is **density < 0.059905026**.

## Checks completed

- Exact reconstruction from the saved degree-63 numerical checkpoint, using closed normalized Laguerre coefficients rather than the producer's recurrence. The explicit rationalization and multiplier guard are checked exactly.
- Independent monomial Fourier identity and positive Fourier normalization.
- Exact harmonic expansion, nonnegative coefficients, angular constant coefficient 1 by the spherical moment formula, and the rank-one origin extension determinant zero.
- Primitive rational Sturm sequence for H on [0,infinity): endpoint/infinity variations 31/31. Physical P on [12,infinity): 30/30. Neither test uses complex roots or floating-point signs.
- Complete pair interval [pi*(13/8),12]: 78 certified subintervals, depth 13. The independent ratio gap and endpoint-secant curvature bound differ from the producer's midpoint-Taylor residual proof. Coverage is checked by exact adjacent endpoints. The rest of the ray follows from P<0 and 0<h<=1.
- Exact Bernstein cap proof: 25 intervals. Exact rational Sturm proves precisely one transition in (1/2,1).
- Full cyclic minimum-side domain: 16934 accepted rectangles, 999 exact sorting exclusions, maximum binary depth 28. Excluded rectangles satisfy 2*t_min*u_min>1. All accepted boxes use rigorous direct or AD mean-value enclosures. Partition area equals 9/2 exactly. There are no unexplained omissions.
- All eight angle/exponent derivatives checked symbolically. The denominator lower bound w^2>=1 is exact. The infinite u>=4 tail and all-dilation theorem are checked as described in `THEOREM_REVIEW.md`.
- Rational Machin pi bounds and squared ball volume give a fully rational strict density rounding. The odd-dimensional volume formula is used, not copied from dimension 15.

The numerical density enclosure is `[0.0599050256891009539149594368291375879349619642325212571560660761301533276111816473100195243573695047372302365585809656105956815021877926547951455272837508 +/- 6.22e-155]`. The objective is evaluated from exact coefficients, not from the optimizer's reported objective. No external kissing or positive-width neighbor-count result is assumed.

## Limits

This establishes one explicit packing upper-bound witness, not optimality or saturation of the method. The comparison is with the latest pinned paper verified in this audit, not a universal literature-completeness theorem. The arithmetic backend is shared python-flint/Arb; the polynomial sign strategy and geometric implementation are separate from the producer. These are internal independent checks, not external peer review.
