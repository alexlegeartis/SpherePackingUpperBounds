# Independent six-dimensional audit

Result: PASS. Strict density upper bound 0.410189987 for unrestricted six-dimensional sphere packings, using the four published average-neighbor bounds as external mathematical inputs.

Frozen discovery input SHA256: 797128e7ada65976ce25926ee84f626a9607208013bfb420243215df2f871606.

The separate verifier `verify_algebra_guards.py` reconstructs Laguerre coefficients using explicit binomial formulas, independently checks the Fourier polynomial, 66 positive rational LDL pivots, exact saturated-block factorization, count cost, mixed-Gaussian Fourier normalization, physical sign repairs, exact Sturm root-free complements, and rational outward rounding. `verify_pair.py` independently certifies the complete pair interval with 1,372 closed leaves and explicitly tests the smaller allowance at each half-open count threshold. These logs and JSON outputs are included here.

The resulting density enclosure is centered at 0.410189986831198867740760436299428501751517129871361804088376415190839315331435110827884942755393531824506321066346, radius below 4.05e-115.

The cyclic and far algorithms were separately source-reviewed, including polynomial identities, positive denominator, exact Fraction domain coordinates, outward binary64 Bernstein subdivision, symmetry exclusions and strict geometric exclusions. The other independent reviewer freshly replayed the complete frozen source: 3,097 far boxes with 414 geometric exclusions, and 153 compact cyclic boxes with 7 symmetry exclusions. These exclusions are legitimate proof steps, not missing domains. This audit is not a clean-room implementation of every continuum checker.

The portable `verification_package/verify.py` was reviewed for bindings to the frozen witness, all exact dimension/count parameters, Gaussian helper, Schur factorization and fresh complete domain subprocess checks. No mathematical blocker was found.

## External dependencies and comparison

De Laat, de Oliveira Filho and Vallentin, https://arxiv.org/html/1206.2608, Section 1.4 and Table 2, define epsilon-tangency edges on [1,1+epsilon) and give six-dimensional average-degree bounds 79,80,81,82 at epsilons .002006,.004024,.006054,.008097. The certificate uses each radius reduced by 1e-6. The published theorem supplies these count bounds; their original SDP witnesses were not rerun. Average degrees enter the ordered-pair sum once per center, so there is no missing factor two.

The pinned paper-reported comparison upper bound is .4102905 in https://arxiv.org/html/2512.04829v2, Table 1. Its source witnesses were not replayed and this audit does not establish an exhaustive worldwide priority claim. The known E6 lower density is pi^3/(48 sqrt(3)).
