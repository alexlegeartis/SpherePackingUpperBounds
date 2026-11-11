# Six-dimensional 99/100 mixture audit

PASS: strict unrestricted sphere-packing density upper bound 0.409135399.

Frozen discovery source: research/density_low_structural_20261111/mix_continue_99hundredths.json, SHA256 41c60fae103355323abdfab8740c64b1089732bc0e0903fd27da845ed4259810. Exact proof inputs are copied alongside this note and verified directly; no assertion depends on convex mixing preserving an optimizer's numerical feasibility.

Distinct arithmetic reconstruction checks explicit Laguerre coefficients, full Fourier identity, all 66 positive rational PSD pivots, exact origin Schur congruence, positive Fourier normalization, exact published-count cost and strict rational outward rounding. Both Gaussian guard weights are zero. Exact Sturm computations prove the full Fourier polynomial positive on [0,infinity) and the physical polynomial negative beyond squared support radius times pi. No sign-repair intervals or second Gaussian width are used in this final witness.

The separately implemented pair Taylor verifier passes 1,672 closed leaves, maximum depth14. All four half-open count thresholds are explicitly checked with the smaller allowance. Fresh full-domain replay of the source-reviewed exact-initialization/outward-subdivision checkers passes 4,359 far boxes (859 geometric exclusions) and 87 compact cyclic boxes (4 permutation-symmetry exclusions). These exclusions are valid proof steps, not omitted domains.

The cyclic summation theorem, origin-coordinate PSD construction, closed support boundary, periodic approximation and external neighbor-count dependencies are identical to the prior reviewed six-dimensional package. Count bounds79,80,81,82 come from de Laat–de Oliveira Filho–Vallentin arXiv1206.2608 Table2 with radii each reduced by1e-6; their original spherical SDP certificates remain external dependencies. This audit is not a formal proof-assistant verification or a clean-room reimplementation of every continuum checker.

The density enclosure is centered at0.4091353980508900600594412074029310240928937752092733306302830174062684616768589755434278958556042087570472399115958 with radius below7.56e-116. Against paper-reported .4102905 (arXiv2512.04829v2 Table1), the outward rounded bound reduces the upper by .2815324751608921% and closes 3.093223388466630% of the additive gap to E6 density pi^3/(48 sqrt(3)). This is a pinned comparison, not an exhaustive priority claim.
