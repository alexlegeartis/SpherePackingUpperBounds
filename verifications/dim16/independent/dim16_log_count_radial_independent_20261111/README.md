# Independent audit of the repaired 16D radial witness

**PASS.** The complete radial gate, combined with the separately audited spherical count and log-hinge geometry, gives

    density16 < 0.024920988.

The exact expression begins 0.0249209879146585407778922410323795054. This is a small improvement over the displayed published comparator 0.02492121, not a substantial gap closure.

## Replay

Use Python 3.12 with python-flint 0.9.0:

    python verify.py
    python negative_controls.py

There are no producer-code imports. Generated reports and full exact Sturm chains are rewritten locally. Assertions are required; both scripts reject `-O`.

## Independent algebra

The degree-63 source is reconstructed from the exact binary64 coefficients using the closed Laguerre coefficient formula, rather than the producer recurrence. A second identity transforms every physical monomial directly: in dimension 16, the Gaussian Fourier map sends X^k to k! L_k^7(Y). The resulting polynomial agrees exactly with the stored H.

The repair adds epsilon(1-X/4), epsilon=1/1000000, to the normalized physical polynomial and subtracts that same degree-one polynomial from H. The subsequent common division by 1-epsilon is applied to P, H, C and the count multiplier. H(0)=1 exactly. All scalars have the required nonnegative signs.

## Complete signs and pair domain

Independent rational Sturm chains have degrees 63,62,...,0. Endpoint and positive-infinity variations are 31 and 31 for H on [0,infinity), and likewise 31 and 31 for P on [21/2,infinity). The endpoint signs are strictly positive and strictly negative, respectively. This does not rely on the producer's interval complex-root isolation.

The finite pair domain is v in [1,7/4], X=2 pi v. A separate first-derivative mean-value proof passes on 59 intervals, compared with the producer's 36 second-order Taylor intervals. Every profile kink is an exact partition endpoint. The count threshold v=5000/4997 is checked from the left with the allowance and from the right without it, including equality. There are no gaps. For v>=7/4, X>=7pi/2>21/2, so the negative physical ray and the nonnegative kernel allowance prove the entire tail, including the later profile kink at 2.

## Rounding and scope

The count input is the separately proved pointwise integer cap 7802. The objective is P(0)+C+7802 lambda; the ball volume is pi^8/(256*8!). Rational alternating-series bounds for Machin's identity enclose pi and prove the strict displayed density bound, independently of Arb's pi routine. ROUNDING.json retains the exact rational endpoints.

Six negative controls check a missing pair interval, wrong count-side semantics, wrong Fourier repair sign, the unrepaired source's rigorously positive pair point, an overstated decimal bound, and disabled assertions. The original source failure remains explicit; it is not hidden by the repaired result.

The complete spherical and profile hypotheses must also be replayed. The aggregate review and theorem linkage are in the sibling dim16_log_count_combined_independent_20261111 directory.
