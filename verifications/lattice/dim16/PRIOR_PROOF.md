# 16D lattice-packing density upper certificate

**Claim:** every lattice packing of equal balls in R^16 has packing density at most **0.018743418752**.

This package combines Zubrilina's modulo-two class bound with a radius-sensitive weighted-frame inequality. It concerns lattice packings only. It does not establish an unrestricted sphere-packing bound or optimality.

Run from any directory:

    python verify_prior.py

Requires Python 3.10+ and SymPy (tested with 1.14). Verification uses rational arithmetic and exact real-root counts. Floating-point optimization is not part of the proof. Runtime depends on the symbolic-arithmetic environment. `prior_certificate.json` is self-contained and the verifier reconstructs the two polynomials from the normalized Laguerre coefficients before checking them.

## Comparison

The safe decimal is about 0.01383% below the previously reported refined bound 0.018746011555. This is a modest improvement, not a substantial reduction of the gap to the best known lattice packing. Worldwide priority must be checked separately. The unrounded rational certified value is recorded in the JSON (approximately 0.01874341875142949).

## Geometric lemma

Normalize the shortest lattice vector to have length 1. Take distinct opposite pairs in a nonzero class of Lambda/2Lambda, represented by v_i, with squared norms t_i in [2,9/4]. Then

    sum_i (10−3 t_i)/4 <= 16.

Proof: if i≠j, both (v_i±v_j)/2 are nonzero lattice vectors. Hence

    |<v_i,v_j>| <= (t_i+t_j−4)/2.

Put a_i=t_i−2, S=sum a_i, Q=sum a_i², and let N be the number of pairs. The Gram rank is at most 16. Its trace and squared trace give

    (2N+S)²/16 <= sum t_i² +(N/2−1)Q + S²/2.

Since 0<=a_i<=1/4, Q<=S/4. Therefore

    C_N(S) = 4N−N²/4 +(4−N/8)S +7S²/16 >= 0.

The ordinary relative bound, applied to normalized representatives with absolute inner products at most 1/9, gives N<=floor(1280/65)=19. For N<=16 the desired inequality is immediate. For N=17,18,19, C_N is strictly increasing on S>=0, while its values at S=4(N−16)/3 are respectively −35/36, −11/9, −3/4. Thus S>4(N−16)/3 and N−3S/4<16. The verifier checks these rational constants.

## Auxiliary function and all class branches

Let p and q be the rational polynomials reconstructed from the normalized generalized Laguerre coefficients d_k:

    p(x)=sum_k d_k L_k^7(x)/L_k^7(0),
    q(x)=sum_k (−1)^k d_k L_k^7(x)/L_k^7(0).

Use h(v)=p(16|v|²) exp(−8|v|²). Fourier scaling yields

    hhat(0)=(2π/16)^8 q(0).

The verifier proves q>0 on [0,infinity), p<0 on [36,infinity), h<=G on every bin between radii 1 and 3/2, and

    h(v)<=M(10−3|v|²)/4,  2<=|v|²<=9/4.

A class without a vector shorter than sqrt(2) therefore contributes at most 16M in opposite-pair units. A class with its unique shorter pair in bin j is bounded by Zubrilina's conditioned formula

    κ_j = g_j + sum_{b>=J}(g_b−g_(b+1)) A(e_(j+1),e_(b+1)).

The certificate uses κ=max(16M,κ_0,...,κ_(J−1)). It does not replace the conditioned branches by an invalid per-class global cap. The zero class 2Lambda has all nonzero lengths at least 2, where h is negative. Poisson summation now gives

    density <= (16/8)^8/8! * [p(0)+2(2^16−1)κ]/q(0).

This final expression is rational and is checked against the safe decimal.

## Exact interval majorants

For each squared-radius bin, a=16 E_b and x−a>=0. The verifier checks

    Gamma_b[1+(x−a)/2+(x−a)²/8]−p(x)>=0.

Since the bracket is at most exp((x−a)/2), this gives h<=Gamma_b exp(−a/2). Heights are rounded upward using Gamma_b divided by the 128th partial sum for exp(a/2), with backward maxima on the long bins. The weighted envelope is checked similarly using

    M [sum_{k=0}^{128}(a/2)^k/k!] (5/2−3x/64)
      [1+(x−a)/2+(x−a)²/8]−p(x)>=0.

All interval checks use exact Bernstein coefficients and subdivision. Global Fourier and tail signs use exact root counts. No sampled positivity verdict is accepted.

## Sources and provenance

- Nina Zubrilina, *Improved lattice-packing bounds in small dimensions*, arXiv:2610.35830, https://arxiv.org/abs/2610.35830 . We retain its short-class caps, Fourier normalization, and exponential-majorant method.
- Original 16D certificate: https://github.com/AstroNinja389/lattice-packing-bound/blob/main/certificates/n016.json ; SHA256 1d6f40a3ca0812758a696750a378dbbafb4dea3f97bdf43f114502aefa43f5b1.
- New degree-61 coefficients were found by a sampled LP with successive continuous-extremum cuts, then independently rationalized and checked. Initially sampled candidates failed Fourier or spatial-tail positivity and were rejected. The final proof depends only on this rational certificate and the geometric lemma.
