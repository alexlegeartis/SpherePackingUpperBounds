# Independent36D Fourier/rank-two/profile audit

PASS: lattice packing density <=0.000027402464246. Exact density is approximately0.00002740246424596642, recorded as a rational in certificate.json. The improvement over the prior safe0.000027428605 is0.0953047157651%. This is lattice-only; no unrestricted or worldwide-priority claim is asserted.

## Frozen inputs and replays

Producer verifier SHA256222678e39357e2c7fea063d69640d052a00fb548368795167b65a4ae5caee65c; candidate6696a83db0890fc5d8faea493c0f5c3fe8bfc05d3cea6ea699d6ce77ccfa3c25; basis323a91456b78dae2632cedc45481d76e7751c7bcc950bce5c6a9d06fc63e7568; certificate27ef83e477b42c420afa0499acd8cb830cfb1f6cd9fd68777b98c8975e89e7a6. All hashes independently matched.

New verifier was independently replayed at1024-bit Arb precision, with an added assertion that exactly three multiplier basis polynomials are used and an exact rational comparison to the stated safe decimal. All240 radial bins/248 cells,193 class profiles, global Fourier and spatial-tail signs, PSD conditions, and rounding passed. The entire certificate JSON equals the producer768-bit output. The separate exact degree83 Laguerre inversion check passes, including unchanged Fourier zero.

The original36D complete certificate was independently replayed in old_audit36_independent/. Its exact polynomial reconstruction,240 radial intervals, four angular-polynomial continuous-domain/code-cap proofs,193 class profiles, even-dimensional normalization, and global signs all passed. This dependency replay is complete, not just a review of producer logs. Independent log: old_audit36_independent/independent_replay.log.

## Angular kernel and Fourier component

The frozen Q37 kernel is reexpanded exactly in the36D normalized Gegenbauer basis (parameter17). All coefficients are nonnegative and the constant is positive. The supplied q equals Q(1)/constant. Q37 hasha57da6cbe33f3f7395ef40b76d3450388115423d60f83270920b10fc4f97212d is checked. Its universal rank-two sign proof is inherited from the completed independent100-decimal degree24 audit, not a new unchecked numerical sign screen.

Dimension36, tau38, normalized Laguerre parameter18 and volume factor(tau/8)^18/18! are consistent. The quadratic Fourier multiplier is repaired to w>=v^2/4+1e-16; this repaired value is included in every recomputed physical/radial check. The differential operator basis and separate exact inverse transform agree. Original Fourier factor is globally positive and the multiplier is positive definite. New physical tail is globally nonpositive.

The radial PSD kernel origin block uses cross z decreased1e-10 and origin z^2/aZ+1e-20. All PSD inequalities pass exactly. The full-line correction uses the already-audited rho=10^-70 geometric-series factor, proved by168^168/168!>10^70. Combined kernel tail coefficients satisfy A>0, B<0 and A+B<0. No fixed-origin transfer from the old36D certificate is claimed: this is a new directly verified joint kernel/Fourier/profile certificate, and the old certificate had no such origin block.

## Reoptimized angular-cut weights and profile budget

The four angular cuts are byte-for-byte equal as JSON data to the original36D cuts. No new angular theorem is assumed. Their Gegenbauer PSD expansions, projective sign domains, cap majorants and spherical-code caps were checked by the completed original independent replay.

The new nu coefficients are nonnegative; hence the old four global cut inequalities may be combined with these new weights. The new long-bin radial heights g are nonnegative and nonincreasing. For cumulative caps C_k made nondecreasing by taking minima over later bounds, the maximum weighted long-bin sum is bounded by sum_k(C_k-C_(k-1))*g_k, by summation by parts. This is exactly the implemented budget. A branch with short-bin j adds g_j minus sum_i nu_i*w_ij; the no-short branch adds neither. Taking the maximum over all193 branches gives kappa. Restoring the global angular-cut budgets adds exactly2*sum_i nu_i*constant_i, outside the factor2*(2^36-1) multiplying kappa. This matches the original cut normalization, and the new objective includes every such term divided by qzero. Negative individual cut weights cause no issue because the global inequalities, not a sign assumption on those weights, are used.

All arithmetic in this class budget is exact rational. The final density uses the complete new physical value at zero, kernel origin and corrected profile budget. Files verify.py, replay.log, certificate.json and check_fourier_transform.py preserve the independently executed proof.
