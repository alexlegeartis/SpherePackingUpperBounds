# Applicability below 17D and above 21D

The periodic cyclic proof, solid-harmonic kernel, and origin Schur reduction are not special to dimensions 17–20. The implementation accepts 3≤n≤48, with the stated finite-support and exact-verification conditions. Applicability does not imply that the small chosen function space beats the literature.

## Valid but noncompetitive examples

The following two examples now pass the complete producer verifier. Replay them, including exact Schur factorization, with:

```sh
python run_range_applicability.py
```

| Dimension | Fully checked raw density enclosure | With exact origin Schur reduction | Published comparator | Outcome |
|---|---:|---:|---:|---|
|16|0.02519483655593481219…|<0.025150156|0.02492121|Valid but about 0.919% weaker|
|22|0.00396025920867277881…|<0.003908297|0.0037259419689206|Valid but about 4.894% weaker|

16D uses seven radial nodes, angular degree 16 and Gaussian-polynomial degree 31. The complete checks cover the pair interval, 56 compact cyclic cells and 21 far cells. Its source polynomial has two forbidden physical-sign islands; the slow Gaussian correction has c=157/25, eta=53/10^16, zeta=0. Exact Sturm complements and Bernstein/Arb bounds repair the intervals [59247/1000,59299/1000] and [14849/100,149581/1000] in X=pi|x|^2. The Fourier-side cost and changed normalization are included. Source: candidate16_wide_guarded.json; verifier: verify_guarded.py, reciprocal_guards.py, slow_gaussian.py; log: audit16_wide_guarded_full.log.

22D uses five radial nodes, angular degree 16 and Gaussian-polynomial degree 39. Exact PSD, Fourier identity, both global one-variable signs, the pair domain, all 20 compact cells and all 10 far cells pass. Some far cells require legitimate impossible-triangle exclusions. No Gaussian guard is needed: both repair coefficients are exactly zero. Source: witness22_degree39/dim22.json; verifier: verify.py; log: audit22_degree39_full.log. The older degree 31 cyclic candidate likewise passed the algebra, global signs and pair checks, but is numerically weaker.

Both examples also passed independent full exact-rational audits, including the Gaussian correction where needed, PSD/Fourier/pair checks, complete domains, Schur factorization and rigorous decimal rounding. The 16D audit proves 77 domains with five exact impossible-triangle exclusions; the 22D audit proves 30 domains with 50 such exclusions. See ../dim16_guarded_independent_20261111/REVIEW.md and ../dim22_cyclic_independent_20261111/REVIEW.md. These are independent computer-assisted checks, not external peer review or formal proof-assistant verification. No new best bound is claimed.

## Benchmarks and coverage

For 16D, use Table 1 of [Tutunov et al., arXiv:2512.04829v2](https://arxiv.org/html/2512.04829v2), which reports 0.02492121. The older 0.0249944093845237 in the indexed Cohn table is not the newest comparison. The same paper reports 15D 0.03372585. We have not independently reconstructed these literature witnesses.

For 22D, 23D and 25D, the pinned comparisons are 0.0037259419689206, 0.0026842798864291, 0.0013841907222857, respectively, in [Cohn's density table](https://cohn.mit.edu/sphere-packing/), citing [Afkhami-Jeddi et al.](https://arxiv.org/abs/2006.02560v2). On 2026-11-11, the search index returned the full table with a four-month-old crawl; direct access returned 502. These are explicit comparators, not an exhaustive current-worldwide-best assertion.

This new cyclic search has explored 15, 16, 17–23 and 25. The 15D numerical result 0.033970452 and 23D 0.002987075 are noncompetitive and not fully certified. The 25D wide trial 0.001532139 is noncompetitive and also fails exact PSD. The 21D trial does not beat the repository's already stronger registered 21D bound. The present series has not tested 6–14 or 26–47. Existing repository improvements in 4D, 5D and 48D come from earlier work and should not be attributed to these new runs. Dimensions 1, 2, 3, 8, 24 are already solved exactly, so a strict improvement over their sharp density is impossible.

A negative result from one finite basis is not a failure of the general method. In particular, changing to solid harmonics also changes each finite radial span; the performance difference cannot be assigned solely to replacing pointwise negativity by cyclic negativity.

Recent lattice-only advances, such as [Zubrilina 2610.35830](https://arxiv.org/abs/2610.35830), do not furnish upper bounds for arbitrary nonlattice packings and are not interchangeable with our general-packing comparators.
