# Six-dimensional patch5 follow-on audit

PASS: strict unrestricted sphere-packing density upper bound 0.409128533.

Frozen discovery source: research/density_low_structural_20261111/patch5_exact.json, SHA256 19bddaa444bb214613549d2ee1ec6c786c585c4bc1c72f7fcd3977acef4ef521. Exact proof inputs are copied alongside this note and verified directly; no assertion depends on convex mixing preserving an optimizer's numerical feasibility.

Distinct arithmetic reconstruction checks explicit Laguerre coefficients, full Fourier identity, all 66 positive rational PSD pivots, exact origin Schur congruence, positive Fourier normalization, exact published-count cost and strict rational outward rounding. Both Gaussian guard weights are zero. Exact Sturm computations prove the full Fourier polynomial positive on [0,infinity) and the physical polynomial negative beyond squared support radius times pi. No sign-repair intervals or second Gaussian width are used in this final witness.

The separately implemented pair Taylor verifier passes 1,866 closed leaves, maximum depth14. All four half-open count thresholds are explicitly checked with the smaller allowance. Fresh full-domain replay of the source-reviewed exact-initialization/outward-subdivision checkers passes 4,457 far boxes (916 geometric exclusions) and 83 compact cyclic boxes (4 permutation-symmetry exclusions). These exclusions are valid proof steps, not omitted domains.

The cyclic summation theorem, origin-coordinate PSD construction, closed support boundary, periodic approximation and external neighbor-count dependencies are identical to the prior reviewed six-dimensional package. Count bounds79,80,81,82 come from de Laat–de Oliveira Filho–Vallentin arXiv1206.2608 Table2 with radii each reduced by1e-6; their original spherical SDP certificates remain external dependencies. This audit is not a formal proof-assistant verification or a clean-room reimplementation of every continuum checker.

The density enclosure is centered at 0.409128532917618682332912747590282332997433318390373354268480707508683656225482605709462523483709564533108876618501 with radius below 2.95e-115. Against paper-reported .4102905, the outward rounded bound decreases the upper by .2832059236077852% and closes 3.111609721596990% of the gap to E6.
