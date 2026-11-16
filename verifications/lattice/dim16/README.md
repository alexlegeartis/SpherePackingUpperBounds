> Reproduction from the repository root: `python verifications/lattice/replay.py --dimension 16`. Large JSON inputs are compressed; the replay checks decompressed hashes in a temporary copy before running the proof programs below.

# Lattice packing in 16D: arithmetic trace-gap certificate

**Safe upper bound: 0.018655117428.** This is for lattice packing density, not arbitrary sphere packings. It is approximately0.471% smaller than the preceding certified0.018743418752. This is a modest improvement, not a percent-level breakthrough or an optimality proof.

Run `python run_all.py` from this directory. Requirements: Python3.10+, NumPy and SymPy. SciPy or an optimizer is NOT needed for verification. The replay checks the preceding bound, finite arithmetic near-frame witnesses (four metric-certificate families and all prime-support branch trees), exact polynomial signs, all spatial intervals, and2498 rational linear-program dual certificates. Runtime depends on the environment. Verification rejects Python's optimized `-O` mode.

## Mathematical structure

Normalize the lattice minimum length to1. The prerequisite certificate proves reciprocal covolume rho<5220<5248. The arithmetic theorem in `NEARFRAME_PROOF.md`, with its exact finite replay, proves this dichotomy:

* rho<=4096; or
* every16 distinct opposite pairs from one nonzero class modulo2Lambda have sum of squared lengths at least32+4/15.

The theorem allows one representative shorter than sqrt(2); this extension is essential. The low-density alternative gives density below0.01476 (using pi<22/7), hence below our claimed bound. We handle the other alternative by auxiliary-function optimization with the full first16-pair trace constraint.

The rational functions are h(v)=p(16|v|²) exp(−8|v|²), with Fourier partner q reconstructed from generalized Laguerre coefficients. The verifier establishes q>0 on[0,infinity), p<0 on[36,infinity), and h<=G on the518 bins in1<=|v|<=3/2. It also establishes the affine envelope

    h(v)<=M(1−(4/5)(|v|²−2)),   2<=|v|²<=9/4.

For an empty class (no shorter-than-sqrt(2) representative), this has opposite-pair budget (16−(4/5)(4/15))M. For N=16 pairs this follows from the trace gap; for N<=15 it is immediate; N=17,18,19 follow from the weighted-frame quadratic C_N in the previous proof. The verifier checks the requisite rational constants explicitly.

## Exhaustive conditioned-class treatment

A class with a short pair in bin j uses its upper squared radius E_(j+1). Let x_b be its number of remaining pairs in long bin b, and put a_b=E_(b+1)−2. The original exact conditioned bounds give

    sum_(b<=k) x_b <= A(E_(j+1),E_(k+1)).

All endpoint arguments here use squared radii. The bins are half-open; polynomial checks include both endpoints and are consequently stronger than needed.

Classes with at most15 total pairs have at most14 long pairs and are bounded by filling the earliest admissible bins. If a class has at least16 pairs, let h be the bin containing its16th-smallest pair, equivalently its15th-smallest long pair. Necessary conditions are

    sum_(b<h) x_b <=15,
    sum_(b<=h) x_b >=15,
    sum_(b<h)(a_h−a_b)x_b <=15a_h+(E_(j+1)−2)−4/15.

The last constraint uses UPPER endpoints. Thus every genuinely feasible configuration satisfies it; lower endpoints would be invalid here. Pairs after the16th may not be used to compensate for a too-short initial16-pair frame.

We relax x_b to nonnegative real numbers. This enlarges the configuration set, so maximizing the weighted contribution remains a safe upper bound. Every possible cutoff bin is enumerated. The2498 stored dual vectors are nonnegative and satisfy C-transpose y>=g exactly. Their rational objective b-transpose y, plus the short-pair height, therefore bounds each profile. No floating-point feasibility or optimizer output is trusted by the verifier. It rejects missing or duplicate cutoff cases.

Only14 short branches need this cutoff treatment; in the others the necessary upper-endpoint trace constraint is redundant or the original caps permit fewer than16 pairs. The verifier recomputes this classification exactly. The common class budget kappa is the maximum of every checked case and the empty-class bound.

## Final density and prerequisites

Poisson summation gives the rational bound

    density <= (16/8)^8/8! * (p(0)+2(2^16−1)kappa)/q(0).

The final value is approximately0.018655117427637372 and is compared exactly to0.018655117428.

`verify_prior.py` and `prior_certificate.json` reproduce the prerequisite density bound0.018743418752. This implies rho<5220 using pi>3.14159; that elementary lower bound follows, for example, from Machin's identity with alternating arctangent bounds. The prior geometric weighted-frame lemma is proved in PRIOR_PROOF.md; its rational constants are also checked in the verifier. Exact alternating-series Machin bounds are checked in verify_nearframe.py.

## Sources

The modulo-two framework, short-class A/B caps, Laguerre representation and polynomial exponential majorants follow Nina Zubrilina, *Improved lattice-packing bounds in small dimensions*, https://arxiv.org/abs/2610.35830 . The new ingredients are the weighted class inequality, arithmetic near-frame exclusion, and first16-pair profile constraints. Comparison/priority checks are separate from certificate validity.
