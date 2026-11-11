# Independent 14D curved-profile audit

Conclusion: the frozen certificate proves density strictly below **0.045003383**.

Run `python verify.py` using Python 3.11+ and `requirements.txt`. Run `python negative_controls.py` for six fail-closed tests. Do not use Python's `-O` or `-OO` modes. No solver, network, producer verifier, or other research folder is required. Generated logs/results may include timing differences; immutable code and mathematical inputs are pinned in `SHA256SUMS`.

The independent implementation uses explicit Laguerre and Gegenbauer formulas, primitive exact rational Sturm for both infinite sign rays, exact Bernstein coefficients for the angular cap, an interval automatic-differentiation cyclic checker, and a pair-ratio secant/curvature proof. This shares python-flint rational/Arb arithmetic with the producer, but does not use its complex-root isolation or geometry implementation. The exact mathematical data are intentionally the same witness. See `AUDIT.md` and `THEOREM_REVIEW.md` for scope and proof.

The pinned published comparator is 0.04501031, from [Table 1 of arXiv:2512.04829v2](https://arxiv.org/html/2512.04829v2#S5). The strict rounded bound improves that comparator by at least 0.015389807%. This is an independent computational audit, not external peer review or a formal proof-assistant verification.
