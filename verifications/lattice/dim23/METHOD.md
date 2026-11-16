# Conditioned nested Gram-moment certificate (23 dimensions)

This certificate has no prior density-bound premise, no quotient classification assumption, and no angular pair-sum cuts. Two full independent audits passed. The portable runner repeats the complete independent verification.

Normalize minimum vector length to1. A nonzero class modulo2Λ contains at most one opposite pair with squared norm below2. Let selected opposite pairs in one class have squared norms2+δ_i. Their representatives satisfy |<v_i,v_j>|≤(δ_i+δ_j)/2 and δ_i+δ_j≥0, because (v_i±v_j)/2 are nonzero lattice vectors.

For N representatives, S=Σδ_i and Q=Σδ_i², Gram rank gives

    (2N+S)² ≤ n [4N+4S+(N/2)Q+S²/2].

For a long-radius cutoff2+d≤9/4, long δ_i lie in[0,d]. If the class has a short representative δ_0=−z with0≤z≤1/96, the other selected δ_i lie in[z,d]. Consequently

    Q ≤ (d+z)S+2z²−(N−2)dz,
    (N−2)z ≤ S ≤ (N−1)d−z.

The projective relative bound gives

    N ≤ 4n(1+d)/(4+4d−(n−1)d²).

This applies with the short representative too: its normalized absolute correlation with a long representative is no larger than d/(2+d). The supplied verifier checks the strict Gram contradiction below S*=(N−n)/α, α=(n−1)d/[4(1+d)], for every integer N>n and every used cutoff, by exact rational polygon vertices and edge convexity. For N≤n, S≥0 when N≥2; the singleton case is immediate. Therefore

    Σ_selected (1−αδ_i) ≤ n.

For an all-long class this has RHS n. When a short pair lies in a half-open bin[a,b)⊂[191/96,2), its guaranteed contribution is at least1−α(b−2); hence the long-pair RHS is n−1+α(b−2). For other short bins use the valid all-long RHS n. Apply this at every long-bin prefix, with upper squared-radius endpoints in the decreasing weights. This is a necessary relaxation, not a lower-endpoint exclusion.

The radial certificate uses h(r)=p(τr²)exp(−τr²/2) and its Laguerre-normalized Fourier polynomial q. Exact coefficient reconstruction and global Möbius-Bernstein signs establish q≥0 everywhere and p≤0 beyond the support. Each radial bin has a rational envelope proved by Bernstein subdivision and a lower Taylor bound for the exponential. Long-bin heights are monotonized upward.

Each of385 possible short-class branches has a nonnegative rational LP dual for the cumulative A/B constraints and the nested weighted constraints. Componentwise domination and its exact objective certify the class budget κ. Poisson summation then gives

    density ≤ (τ/8)^(n/2) / Γ(n/2+1) · [p(0)+2(2^n−1)κ]/q(0).

Odd-dimensional normalization is checked by squaring and a rational lower bound for π obtained from the alternating arctangent series and Machin's identity. The final decimal is rounded upward.

Reproduction: `python run_all.py` after installing SymPy. Public `python-flint` is an optional accelerator. The default uses exact Möbius-Bernstein global Fourier signs when it is installed, otherwise exact SymPy root counts; tail signs use exact Bernstein subdivision plus a positive shifted-power tail. No private paths or source-density premises are needed.

Class bins are half-open [E_j,E_{j+1}); squared norm 2 is long, and the final endpoint belongs to the nonpositive tail. Polynomial envelope checks on closed intervals are conservative. This convention also justifies A=0 when the sum of upper endpoints is at most 4: actual selected norms are strictly below those endpoints.

This builds on Zubrilina's mod-2 class-profile radial framework (https://arxiv.org/abs/2610.35830), adding conditioned nested Gram-moment profile inequalities. It does not use pair-sum angular cuts.
