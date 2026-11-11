# Independent diagnostic22D audit

This is a noncompetitive test of the general cyclic method outside dimensions17–20. It is not a new record and should not replace the published22D upper bound0.0037259419689206.

The frozen source uses dimension22, Gaussian polynomial degree39, five radial hats, and angular degrees0–16. Exact exclusion square is11/4, and exact truncation square33/10. No mixed-Gaussian correction or extra packing assumption is used.

algebra_pair.py independently checks17 PSD blocks with86 positive rational pivots, the direct monomial Fourier identity using parameter n/2−1=10, global Fourier positivity and physical tail via exact Sturm sequences, and the solid-harmonic pair inequality by Arb Taylor enclosures. The latter needs18 leaves of maximum depth6. The density interval is centered at0.00396025920867277881114583224676633441546303387431062912742043002084932925236870904698415416362245938264656048664040. Rational Machin pi bounds certify the outward decimal0.003960260.

domains.py reconstructs normalized solid harmonics by their explicit rising-factorial formula. Exact Bernstein subdivision verifies20 cyclic and10 far cells. It excludes a far box only when exact rational endpoints satisfy w_min>u_max+v_max and(w_min−u_max−v_max)^2>4u_max*v_max, proving the triangle inequality impossible. Equality is never excluded. Floating point is used only to choose the subdivision axis. Full verification passed all30 domains, with50 exact infeasible-triangle exclusions. See full.log and domains_audit.json.

The finite PSD cyclic-sum argument, periodic Poisson bound and periodic approximation are identical to the audited narrow17–20D certificates. Its mathematics is dimension-generic; the Fourier parameter, harmonic parameter and ball volume are evaluated in22D here. This provides a concrete applicability test, not a guarantee that optimization improves known bounds in every dimension.

schur.py independently solves the origin block and verifies the exact PSD factorization after maximal origin subtraction. Exact rational rounding proves a tightened bound0.003908297, still worse than0.0037259419689206. It depends on the full source-domain checks and is not independently sufficient on its own.
