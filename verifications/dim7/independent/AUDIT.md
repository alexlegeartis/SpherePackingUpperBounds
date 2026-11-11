# Independent seven-dimensional density audit

PASS: unrestricted sphere-packing density <0.320584772.

Frozen source: research/density_low_structural_20261111/patch_dim7_exact.json, SHA2561cf76290693c7461ad4da8e5389db283ea5bce02c152ea7ab0115bbf88948bf1. Exact proof inputs are copied alongside this note. The artifact verifies the final exact function directly, not the optimizer's reported feasibility.

The separate algebra verifier constructs degree31 Laguerre polynomials from generalized binomial coefficients for n/2=7/2 and independently verifies the full Fourier monomial transform. All66 rational PSD pivots are positive; the exact origin Schur factorization and recorded objective agree. Gaussian correction weights are zero. Exact Sturm computations establish global Fourier positivity and physical tail negativity. The volume of the radius1/2 ball in dimension7 is pi^3/840; rational Machin bounds certify the strict outward decimal.

The resulting density interval is centered at0.3205847715830073938695391744713529466970446527646514319059435920441989479862831315082043475569821562950856094394234 with radius below6.93e-116.

Independent pair Taylor verification resolves2,727 closed leaves, maximum depth15, and explicitly tests all five count thresholds using the smaller right-hand allowance. Fresh complete replay of the source-reviewed far and compact cyclic algorithms resolves6,487 far boxes with1,223 valid geometric exclusions, and81 cyclic boxes with2 valid symmetry exclusions. Exclusions are justified proof steps, not missing domains. The algorithms construct the dimension-specific normalized Gegenbauer basis from n=7.

The existing periodic PSD/count summation argument, closed kernel support and compact/far partition apply unchanged: p4, K12, R3/2 and minimum separation1. This is not a proof-assistant formalization or a clean-room reimplementation of every continuum checker.

## External inputs and comparison

Live inspection of de Laat, de Oliveira Filho and Vallentin, https://arxiv.org/html/1206.2608 Table2, confirms dimension7 bounds136,137,138,139,140 at epsilons .001001,.002006,.003013,.004024,.005037. The certificate uses each corresponding radius reduced by1e-6 and preserves the half-open neighbor convention. The published SDP bounds remain external mathematical inputs, not independently recomputed certificates.

The paper-reported comparator .32111475 was checked in https://arxiv.org/html/2512.04829v2 Table1. Its witness is not replayed here. E7 lower density is pi^3/105 = .2952978731457125730997744292104894781… . Against that comparator, the outward new bound decreases the upper by .1650431816040839% and closes2.052835449428060% of the additive gap. This is a pinned comparison, not an exhaustive worldwide priority claim.
