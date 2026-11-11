# Independent audit of the 16D log-hinge cyclic geometry

**PASS.** This directory independently verifies the frozen angular/profile component. It does not alone assert a density bound.

Run `python verify.py` using Python 3.12 and python-flint 0.9.0. Run `python negative_controls.py` for four fail-closed checks. Neither script imports a producer implementation. Assertions are mandatory and `-O` is rejected.

The exact dimension-16 Gegenbauer polynomial is independently reconstructed by its closed coefficient formula, with parameter 7, rather than the producer recurrence. Rational Sturm proves the cap margin and the unique sign transition. All fourteen exact printed profile coefficients are nonnegative; both infinite-tail comparisons pass at 256-bit Arb precision.

The complete saved binary partition is reconstructed from the original rectangle, not accepted merely because its areas sum correctly. There are 31,138 accepted rectangles and 327 strict sorted-domain exclusions; the entire exact area is 9/2. Every accepted rectangle passes a fresh analytic-gradient enclosure, with zero additional subdivisions. This uses endpoint-monotone profile ranges and explicit one-sided slope hulls instead of the producer's automatic-differentiation positive-part implementation.

At profile kinks, the function is locally Lipschitz. The slope hull contains every almost-everywhere derivative on the box, so the coordinatewise fundamental theorem of calculus justifies the centered mean-value bound. Values and derivatives are not evaluated using an arbitrarily chosen derivative at a kink.

For all-scale extension, v phi'(v) is nondecreasing, including upward jumps. Hence the two dilation ratios and their quotient decrease. The angular polynomial is negative up to its one root; the sorted remaining angular values cannot have signs (+,-). The mixed (-,+) case is handled by factoring the middle ratio: the remaining bracket decreases, and after it becomes negative the whole cyclic expression is already negative. This is valid for the full signed sum, without replacing negative contributions by their absolute values.

For the u>=4 tail, sortedness ensures the middle-cosine numerator is positive. The final exact linear slope is at least 7; the positive terms are bounded by B(exp(-3 gamma)+exp(-14 gamma))<1/20000. The largest-angle term contributes at most -1/20000. Thus there is no uncovered infinite radial region.

Input hashes and the complete check counts are in VERIFIED.json. The negative controls remove a rectangle, misclassify a valid box as sorted-out, introduce a negative profile coefficient, and invoke Python optimization; all four are rejected.
