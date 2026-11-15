# Kernel-link review

Review scope: exact identity of the stage 1 assembled exterior kernel and the
mathematical formula whose lower envelopes stage 5 certifies, on the full
`B_1` and `B_5` domains. This is not a rerun of the global density certificate.

The recovered partial implementation was completed and rerun by its recovery
author on 2026-11-13. This is author review, not an independent code review.
The portable command is:

```
python3 verify/kernel_link/run.py --out NEW_DIRECTORY
```

The final run exited 0 in 111.24 seconds, with peak owned-process RSS
298,975,232 bytes. It checked 25 ordered radial cells, the 15 underlying stage 1
exterior branches, 213,550 exact polynomial coefficients, 1,655 nonpolynomial
monomials, 1,950 cap-denominator identities and 25 endpoint faces. All five
in-memory corruptions were rejected, including a nonintegral polynomial
denominator. The small result record is retained in `REVIEW_RESULT.json`.

The reviewer inspected the sign and factor conventions in `build_signed`, the
exterior coordinate change, the full-pair gauge, and the stage 5 `minusK_lb`
formula. The symbolic coefficient identity has no floating tolerance. Cap,
tent, allocation and indicator functions retain their exact function arguments;
these cancel structurally. The result covers the closed angular endpoint by a
separate identity, and covers radial breaks using the closed-inner convention,
without assuming continuity across a break.

The code that converts the interval-bound helpers into their exact mathematical
functions is an explicit, trusted adapter. Its formula-to-code correspondence
and lower-envelope soundness remain written proof obligations discharged by the
formulas in `README.md` and the stage 5 `kernel_box_method.md`; these adapters
have not been formalized in Lean. The 375 singleton component checks are useful
regressions but are not used as a proof of equality. The C++ triangle sign
subdivision and the stage 5 Arb leaf replay must still pass in the full driver.
