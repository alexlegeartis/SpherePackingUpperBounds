# Independent 5D long-edge lens density audit

PASS: unrestricted sphere-packing density <0.509996147.

Frozen exact source SHA256:235920216afbd10b38be5efd3672e9d94b5ec5af8fc54aca3a63d77d9997721f. The copied budget matches the producer's frozen step_allowances_005.json. The proof uses elementary centroid caps only, not later Rankin or spherical-code variants.

## Analytic counting audit

For a long pair at distance d, every common neighbor within short radius r=27/25 lies in the midpoint ball of squared radius rho²=r²-d²/4. Pairwise center distances are at least1. The centroid identity gives M(M-1)/2 <= M² rho², hence M <=1/(1-2rho²) when rho²<1/2. Integrality gives the recorded conservative caps83,8,5,3,2,1 on the corresponding length ranges. Every individual bin's cap was recomputed rationally from its lower endpoint; no external kissing input enters this lens cap.

The nonnegative step allowance eta depends only on long-base length. An eligible triangle has two sides <=r and a base at least41/25>r, so it has a unique eligible base. Summing triangle allowances is therefore bounded by sum_pairs B(d)eta(d). Local root-PSD expansions count each unordered triangle's cyclic kernel sum twice. The Fourier/pair expansion counts each unordered pair twice as well. Thus the verified physical-tail inequality f(d)+B(d)eta(d)<=0 pays exactly once for the permitted positive far contribution, with no additional origin objective term and no missing factor2 or3.

For periodic packings, enlarge the period lattice beyond the finite kernel/lens interaction diameter. Count physical centers and triangle/pair orbits modulo simultaneous lattice translation; translated copies of motif centers are not omitted. The same finite incidence identity and caps apply. Schwartz Fourier sums converge, the kernel/lens sums are locally finite, and the standard large-patch periodic approximation then gives unrestricted packings. No assumption of lattice packing is introduced.

Outside the lens support the old far sign is retained. On the kernel support boundary and short-radius boundary, closed interval verification also covers the stronger limiting side. At radial step boundaries, eta is the minimum of adjacent allowances; at the two support endpoints it is zero. No pair can receive two overlapping budgets. The ordinary near-neighbor counts have disjoint radial support from these long-edge terms.

## Independent numerical proof

The algebra audit independently reconstructs the n/2=5/2 Laguerre polynomial, direct Fourier monomial transform, all66 positive rational PSD pivots, exact Schur saturation, positive normalization and rational outward rounding using volume pi²/60. Global exact Sturm signs pass. Gaussian corrections are zero.

Independent pair verification resolves1,074 closed leaves, with all four published-count endpoints checked explicitly at the smaller allowance. Fresh compact cyclic replay resolves115 boxes (6 valid symmetry exclusions). The compact checker is source-reviewed and freshly replayed, not independently reimplemented here.

The new budget checker uses derivatives with respect to X=pi*d² rather than the source checker's radius derivatives. It proves all104 CLOSED bins using113 interval Taylor leaves, maximum depth2. It binds the exact radial coefficients, checks all centroid caps and explicitly verifies the minimum-adjacent boundary values. It does not trust stored success flags or numerical construction minima.

The independently written new far checker reconstructs exact tensor Bernstein coefficients using explicit Gegenbauer formulas, performs outward-rounded subdivision, and derives exact squared-distance extrema on each box. It resolves2,879 boxes:915 passing leaves including11 genuinely relaxed leaves, and526 valid geometry exclusions, maximum depth22. Radial initialization covers the three symmetric rectangle types split at27/25. A box receives at most the minimum eta among every closed bin intersecting its exact distance range, and receives zero if it meets outside-support or outside-short-radius regions. Thus discontinuities and both one-sided limits are covered conservatively.

The exact witness at short sides1,1 and long base1883/1000 has positive kernel value approximately4.303908118922817e-7. The old pointwise far-sign condition is genuinely violated; the new theorem, not merely a tighter old search, is required. Its exact rational value is preserved in positive_far_witness.json.

Density center:0.509996146405708227801189161260425467393448829426584838558811653019445113902227639581328304406372639760404759668067, radius below2.24e-115.

## External assumptions and comparison

The near-neighbor counts45,46,47,48 remain published inputs from https://arxiv.org/html/1206.2608 Table2, at thresholds reduced by1e-6. Their spherical SDP certificates are not rerun. No additional external code bound is needed for the lens caps.

Versus paper-reported .5122645 (https://arxiv.org/html/2512.04829v2 Table1), the outward bound reduces the upper by .4428089395224537% and closes4.825575909597056% of the gap to D5 density pi²/(15 sqrt2). The comparator witness is not replayed. This is a computer-assisted audit, not formal proof-assistant certification or a literature-priority claim.
