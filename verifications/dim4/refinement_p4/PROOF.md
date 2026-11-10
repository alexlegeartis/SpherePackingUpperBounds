# Four-dimensional radial-degree-four refinement

The same argument and four external neighbor-count bounds as the [parent proof](../PROOF.md) apply. Replace w(r) by (T_0(z),...,T_4(z)); all other definitions, angular degree 16, cutoff 333/200, and Fourier degree 31 are unchanged. The exact rational matrices and finite-decimal rational coefficients are frozen here.

The resulting unrestricted density is strictly below **0.634526023**, with Arb enclosure centered at 0.6345260220720713322566555674163184830382171356283089733216206914541048187119 and radius below 7e-77. This closes approximately 8.189832% of the pinned gap from 0.63610277 to pi²/16, not 20%.

Producer replay passes: all 17 exact PSD blocks; exact Fourier/tail Sturm; pair 231 boxes, depth 8; far 135399 boxes, depth 36; cyclic 127 boxes, depth 10. The check is global, not sampled feasibility. Run python verify.py without disabling assertions. Independent audit status is recorded separately. The same IEEE-754/NumPy/Arb trust and external neighbor-count assumptions as the parent package apply. No external peer-review, optimality, or priority claim is made.
