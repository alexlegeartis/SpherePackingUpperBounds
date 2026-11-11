# Independent five-dimensional guarded density audit

PASS: unrestricted sphere-packing density <0.510129239.

Frozen source: research/density_low_structural_20261111/patch_dim5_guarded_1e8.json, SHA25679881d16bc8d837d13222c0d3c5ca6d47264e031e09cd46738275198d4da85b1. The unguarded predecessor patch_dim5_exact.json has an invalid far-sign lobe and is not certified by this report.

The exact final proof inputs are copied here. Independent generalized-binomial Laguerre construction for n/2=5/2 verifies the full Fourier monomial identity. All66 rational PSD pivots, exact saturated-origin congruence, count objective and normalization pass. Both Gaussian repair weights are zero. Exact Sturm checks prove global Fourier positivity and physical tail negativity. Rational Machin bounds and the volume factor pi^2/60 establish the strict outward decimal.

Density center:0.510129238656660990564535665135547517973647651982139698896619983732308355495665512108885511577228162237453462871390, interval radius below3.98e-115.

Independent complete pair verification passes1,085 leaves with all four strict count thresholds explicitly tested at the smaller allowance. Fresh complete replay of source-reviewed domain code passes2,187 far boxes (348 valid geometric exclusions) and57 cyclic boxes (6 symmetry exclusions). These exclusions are justified proof steps, not missing regions.

## Kernel correction

The corrected kernel is verified directly above; its construction is additionally audited. The dimension5 guard coefficients are reconstructed using an independent explicit Gegenbauer formula. Exact Sturm arithmetic proves F(t)<-1/2 on[-1,7/9]. Its origin block has nonnegative diagonal and determinant zero, and its higher coefficients are nonnegative. With D=F(1), cross=-D/2, the pair correction2cross+D is exactly zero.

An independent rational matrix-entry comparison confirms that the frozen corrected source is exactly patch_dim5_exact plus1e-8 times this guard: only the prescribed origin/constant-radial entries change. Radial and count coefficients are unchanged. Final origin saturation is recomputed from the corrected block, so the guard's cost is included. The angular cap7/9 covers every surviving separated pair on the closed annulus[1,3/2]. No unsampled sign assumption is used.

## Inputs and comparison

Published counts45,46,47,48 at epsilons .003013,.008097,.013259,.017446 are from de Laat–de Oliveira Filho–Vallentin, https://arxiv.org/html/1206.2608 Table2. Each certificate radius is1e-6 below the corresponding published threshold. The published spherical SDP certificates remain external mathematical dependencies.

Pinned paper-reported comparator .5122645 appears in https://arxiv.org/html/2512.04829v2 Table1; its witnesses are not replayed here. D5 lower density is pi^2/(15 sqrt(2)). The outward new bound reduces this pinned upper by .4168278301541489% and closes4.542442927666954% of the additive gap. This is not an exhaustive priority claim or proof-assistant formalization.
