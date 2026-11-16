# Near-frame arithmetic exclusion, 16D

## Claim

Let Λ have minimum length1 and reciprocal covolume ρ. Using the already verified upper bound ρ<5220, we prove:

**Either ρ<=4096, or every16 distinct opposite pairs in one class modulo2Λ with squared lengths t_i>=2 obey sum_i(t_i−2)>=16/105.**

This is a structural constraint, NOT by itself a new universal density bound. The low-density alternative equals the Barnes–Wall16 packing density, approximately0.0147081643974. No global optimality assertion is made.

## Frame metric and determinant

Let V have these16 representatives as columns, G=VᵀV=2I+E, a_i=t_i−2>=0, S=sum a_i, Q=sum a_i². The midpoint lattice condition implies |E_ij|<=(a_i+a_j)/2. Consequently

* tr E=S and tr E²<=8Q+S²/2<=(17/2)S².
* λ_max(E)<=(S+4sqrt(Q))/2<=5S/2. Indeed |E| is entrywise dominated by (a1ᵀ+1aᵀ)/2, whose top eigenvalue is the stated value.
* For S<=1/4, ||E||<=3S and log det(G/2)>=S/2−(17/10)S²>=3S/40>=0. This follows from log(1+x)>=x−x²/[2(1−r)] for |x|<=r, with r=3/8. Thus det G>=2^16 and V is nonsingular.

The lattice L0=V D16/2 is contained in Λ, has covolume sqrt(det G)/2^15>=1/128, and its index I=[Λ:L0] is an integer. Hadamard/AM–GM gives

    I <= (ρ/128)(1+S/32)^8.

If S<16/105, the known ρ<5220 makes I<43. Only primes<=41 can divide it. The metric distortion from D16/sqrt2 to L0 is at most1+5S/4.

## Finite torsion witnesses

For a vector a/p modulo D16 (a∈D16, p odd), center each coordinate residue as b_i/p. Let R²=sum b_i²/p² and M=max|b_i|/p. The squared distance to D16 is

    R² + (sum b_i mod2)*(1−2M).

This follows by independently rounding to Z16 and, if the resulting integer vector has odd coordinate sum, changing the nearest coordinate by1. Because p is odd and sum a_i is even, the parity needing repair equals sum b_i mod2.

Signs of the centered residues do not affect this expression. For p=5,7,9,11, the verifier enumerates all histograms of absolute residue types with total multiplicity<=16. For p=9 it excludes histograms without a unit coordinate, which would represent order at most3. The maximum over histograms of the minimum over nonzero multiples gives unscaled squared norms

    p5:42/25, p7:74/49, p9:122/81, p11:16/11.

For p13,17,19,23,29,31,37,41, the verifier bounds the average over multiples by enumerating only the parity support of the histogram. If E is the set of residue types with odd multiplicity, the number m of nonzero coordinates is at most16 when |E| is even and at most15 otherwise. Its averaged square contribution is m(p+1)/(12p). The repair parity depends only on E, while its maximum absolute coordinate is at least the maximum among E. Extra even multiplicities can only reduce the repair cost. These observations give the exact integer numerator implemented in the verifier. Empty E has no repair cost. Every resulting bound is at most42/25.

After scaling D16 by1/sqrt2, all these nonzero torsion cosets therefore contain a witness of squared length at most21/25. Metric distortion preserves length<1 whenever S<16/105, since

    (1+5(16/105)/4)*(21/25)=1.

Thus I has no prime factor>=5 and no element of order9.

If its3-torsion had rank at least2, choose a2-dimensional subgroup. The map D16/3D16→F3^16 is an isomorphism (3 is odd). The eight nonzero words in this subgroup have total Hamming weight at most16*6=96, so one has weight at most12. Its unscaled squared D16 distance is at most14/9, also below42/25. This excludes3-torsion rank2 under the same distortion bound.

Therefore I=2^a or3*2^a. Since I<43, I<=32. Together with det L0>=1/128, this gives ρ=I/det L0<=4096, the claimed dichotomy.

## Replay

Run `python verify_nearframe.py` from workspace root. It requires NumPy and uses integer arrays and exact rational comparisons, not floating-point feasibility. `nearframe_audit.json` records each finite bound and the rational thresholds. The geometric arguments above are not replaced by the numerical enumeration.

## Audited extension: a short representative is allowed

The restriction t_i>=2 can be removed. The midpoint condition gives a_i+a_j>=0 for each distinct pair. Hence at most one a_i is negative. If it equals−z, the other15 are at leastz, so S>=14z>=0. For fixed S,z the maximal sum of squares is

    Q <= (S−13z)²+15z² <= S²,  0<=z<=S/14.

The Frobenius bound therefore remains unchanged. For the top eigenvalue bound, one must NOT claim entrywise |E|<=C on the diagonal, where C=(a1ᵀ+1aᵀ)/2. Instead, for any real x,

    xᵀEx <= |x|ᵀC|x| <= λ_max(C)||x||²,

because E and C have equal diagonal entries and the off-diagonal C_ij are nonnegative. Thus λ_max(E)<=(S+4sqrt Q)/2<=5S/2 still holds. The determinant and torsion arguments follow unchanged.

Consequently the proved dichotomy concerns **any16 distinct opposite pairs in one mod-2 class**, including a shorter-than-sqrt2 pair:

    either ρ<=4096, or sum_i |v_i|² >=32+16/105.

This extension was independently cross-audited with the dual and structural workers.
