# 18D near-frame exclusion with admissible order5 torsion

**Status:** exact finite checks and independent analytic/coordinate audits passed. This is a conditional geometric constraint, not by itself a new universal packing bound.

## Statement

Let Λ⊂R^18 have minimum length1 and reciprocal covolumeρ<33100. The published lattice density upper boundΔ18≤0.0103704834 implies this premise. Let v1,…,v18 represent18 distinct opposite pairs from one class modulo2Λ, and put S=Σ(|vi|²−2).

Then eitherρ≤32768 or S>1/25. Individual representatives belowsqrt2 are permitted.

## Metric and index

Setδi=|vi|²−2 and T=[v1 … v18]/sqrt2. Midpoint constraints giveδi+δj≥0 and |(TᵀT)ij|≤(δi+δj)/4. These imply ||δ||₂≤S, including the possible single negativeδi, and

    ||T||²≤1+(1+sqrt18)S/4≤1+(21/16)S.

Assume S≤1/25. The dimension-uniform determinant lemma (`DETERMINANT_PROOF.md`) applies because S<8/17, provingdetT≥1. Hadamard givesdetT≤(1+S/36)^9. Since T(D18/sqrt2)⊂Λ, its index I is a positive integer with

    I=ρdetT/256 <(33100/256)(901/900)^9<131.

The metric distortion is at mostD=421/400. WriteΓ=sqrt2 T⁻¹Λ⊃D18. Any nonzero element ofΓ with a representative whose physical squared length in D18/sqrt2 is w and Dw<1 is forbidden. We classify the possible odd part ofΓ/D18.

## Residue formula

An element a/p moduloDn, with a∈Dn and p odd, has centered coordinate residues bi/p. Its squared distance toDn is

    Σbi²/p² + (Σbi mod2)(1−2max|bi|/p).

All quantities depend only on absolute residue multiplicities, not their signs. Physical squared distance forDn/sqrt2 is half this expression. The formula follows by rounding toZn, then repairing the coordinate-sum parity if necessary. It applies to all nonzero multiples of the element.

## Small torsion

Complete absolute-residue histogram enumeration, including zero coordinates and every nonzero multiplier, gives:

* order7: physical witness≤43/49;
* order9: physical witness≤22/27;
* order11: physical witness≤95/121.

For order9 a unit coordinate is required, which exactly excludes elements of order at most3. These witness bounds all remain<1 after multiplication byD.

A rank-two3-torsion subgroup corresponds to a2-dimensional ternary subspace ofF3^18. Its8 nonzero words have total Hamming weight≤18·6=108. Thus one has weight≤13. A ternary vector of weight m has physical D18/sqrt2 squared distance m/18 when m is even, and(m+3)/18 when m is odd. For m≤13 this is≤8/9, again forbidden underD. Therefore the3-primary part has order at most3.

## All other possible primes, in quadratic rather than exponential time

Only primes≤127 can divide I<131. Consider a prime p≥13 and the parity supportE of the multiplicities of the absolute residue types1,…,(p−1)/2. Writeh=(p−1)/2.

- If E is empty, no parity repair is needed, and the mean unscaled squared distance over nonzero multipliers is at most18(p+1)/(12p).
- If |E|=1, at most17 coordinates are nonzero. The repair parity depends on that sole residue type. Bounding the maximum absolute coordinate below by that type gives an average correction at most(p+1)/(4p). Thus the mean is at most17(p+1)/(12p)+(p+1)/(4p).
- If |E|≥2, select two distinct residue types. Normalize their ratio to1,a with2≤a≤h. Dropping the parity indicator only increases the nonnegative repair cost, so its mean is at most

      1−2·average_k max(|k|_p,|ak|_p)/p.

  Add18(p+1)/(12p) for the squared-coordinate contribution.

For each primep13,…,127, the verifier checks all a and k exactly. Every resulting physical bound is≤37/39. This covers every parity support without enumerating2^h supports. Moreover

    (421/400)(37/39)=15577/15600<1.

Thus no primep≥13 divides I. Together with the order7/order11 exclusions, only3 and5 remain among odd primes.

## Retaining order5, but excluding order25 and rank-two5-torsion

In18D, order5 cannot be excluded: the histogram with9 coordinates of absolute residue type1 and9 of type2 has minimum physical squared distance1. Exact enumeration proves it is the only exception: every other histogram has a witness≤22/25. SinceD(22/25)<1, every nonzero order5 element inΓ/D18 must therefore have the9+9 histogram, in particular full Hamming support18 modulo5.

A rank-two5-torsion subgroup is impossible. Its24 nonzero words would all have weight18, so total weight432. Each of18 coordinate functionals is nonzero at at most20 of the24 words, giving total weight≤360, a contradiction.

Suppose an element a/25 has order25. Its order5 multiple a/5 forces all18 coordinates of a to be units modulo25. Averaging over the20 unit multipliers, each centered coordinate has mean squared value21/250. Since18 is even, the parity support of absolute unit-type multiplicities is either empty or has at least two members.

If empty, the repair cost is zero. Otherwise normalize two odd unit types to1,a, where a ranges over2,3,4,6,7,8,9,11,12 up to sign. Exact checking gives

    average_{k unit mod25} max(|k|_25,|ak|_25)/25 ≥42/125.

Consequently the unscaled mean distance is at most

    18·21/250+1−2·42/125=46/25,

or physical squared distance≤23/25. SinceD(23/25)<1, order25 is impossible. Therefore the5-primary part has order at most5.

## Conclusion

The odd part ofI divides15, soI=2^a3^b5^c with b,c∈{0,1}. Among such positive integers below131, the largest is128. Therefore

    ρ=256I/detT≤256·128=32768.

This proves the stated alternative. The finite computations and all rational thresholds are verified by `verify_small_torsion.py` and `verify_remaining_torsion.py`.
