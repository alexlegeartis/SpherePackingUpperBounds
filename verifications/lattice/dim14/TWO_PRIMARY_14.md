# The 14D trace-gap theorem used by this certificate

Assume a lattice Λ has minimum nonzero vector length at least1. Write ρ=1/covol(Λ). Subject to the prior bound Δ14≤.03525, either ρ≤512 or every14 distinct antipodal pairs in one class modulo2Λ satisfy

    Σ|v_i|² >=28+12/25.

The prior input is supplied by the included original radial certificate, whose density claim is .0352392080. `run_theorem.py` checks its dimension and threshold, replays it, then runs the arithmetic dependencies below. The final radial certificate and low/high branch composition are checked by the package's final verification runner.

## Near-frame normalization

Let G=VᵀV, a_i=G_ii−2, and S=Σa_i. Assume S<12/25. Midpoints give a_i+a_j≥0 and |G_ij|≤(a_i+a_j)/2 for i≠j. The determinant lemma in `DETERMINANT_DIMENSION_TRANSFER.md`, checked by `determinant_transfer.py`, gives det(G/2)≥1 since12/25<8/13. It includes the possible single negative a_i.

Also G/2≤D I with D≤1+19S/16<157/100; this follows from sqrt14<15/4 and the signed-excess Rayleigh estimate. Set L0=(V/2)D14⊂Λ and H=Λ/L0. Its index satisfies

    I < (964/64)(1+(12/25)/28)^7 <17.

Here the original density bound, exact Machin lower bound for π, and Hadamard imply ρ<964 and the displayed index bound. All these rational comparisons are replayed by `verify_theorem.py`.

## Canonical coordinates and even-modulus parity

An order-p point is x=a/p modulo D14 with a∈D14; the physical vector is Vx/2. Write a=b+pz, retaining the parity lift ε=Σz mod2. After centering kb coordinatewise to b(k), its correction parity is

    ε_k=kε+Σ_i(kb_i−b_i(k))/p mod2.

For even p, this parity must not be inferred from Σb_i(k). The metric verifiers explicitly retain both parity lifts. Coordinate sign changes preserve D14; absolute-residue histograms and both lifts cover every case.

Every actual physical vector gives xᵀGx≥4. The finite metric models use these inequalities, the pairwise Gram bounds above, and a_i≥−1, then average within coordinate groups of equal residue. Positive semidefiniteness is omitted only as a relaxation. Rational duals prove lower bounds on S. The verifiers regenerate every row, all histograms, and all parity cases, then check the dual signs, equalities, and objective values exactly.

## Binary subgroup and the integral kernel

The coordinate-residue map on H[2] has at most the integral odd-sum coset as its kernel. That coset contains each e_i modulo D14. Minimum1 would imply G_ii≥4 for every i, hence S≥28. Therefore the kernel is absent here.

For nonzero binary words, the exact sign-averaging argument in `BINARY_TRACE_14.md` excludes weights2,4,6 when S<4. Consequently all nonzero binary words have weight≥8. Rank4 is impossible: its15 nonzero words would have total weight≥120, while14 coordinates contribute at most8 ones each, total≤112. Thus rank H[2]≤3.

## Order8 and order4

`n14_even8_metric.json` contains2856 primitive order8 histogram/lift cases. Their exact trace lower bounds are at least28/15, so order8 is excluded. An all-even numerator of an exact order8 point would have a nonzero multiple in the already excluded integral odd-sum kernel. Thus the primitive restriction omits no surviving case. Larger2-power orders contain order8 and are excluded as well.

`n14_even4_metric.json` contains96 order4 cases except the two parity lifts of the good histogram

    (zero coordinates, odd coordinates, coordinates2 mod4)=(0,8,6).

Every bad case requires S≥2/3. All-even exact order4 numerators reduce to the integral kernel already excluded. Thus every surviving order4 point has the good histogram.

If x has order4 and y∈H[2], let c be the8-coordinate odd support of x and d the binary image of y. If d has a coordinate outside c, x+y gains a zero coordinate modulo4 and contradicts its required good histogram. Thus support(d)⊂c. Binary minimum8 implies d=0 or d=c. Injectivity gives y=0 or y=2x. Hence an extension containing order4 has binary rank1; with order8 excluded, its2-primary subgroup is cyclic of order4. Otherwise the2-primary subgroup is elementary abelian of rank≤3. In either case its order is≤8.

## Odd primes and conclusion

The exact metric minima for primes3,5,7 are respectively12/25,39/25,31/20, from14,119,679 cases in the three included prime JSON files. They exclude these primes under the strict assumption S<12/25, including equality of the first minimum with the threshold.

For primes11 and13, `canonical_odd14.py` enumerates all residue histograms and gives canonical physical squared witnesses73/121 and8/13. Multiplication by157/100 keeps each strictly below1, excluding these primes. No other odd prime is below17.

Thus I is a2-power and I≤8. Since covol(L0)≥1/64, ρ=I/covol(L0)≤512, proving the asserted alternative.

## Scope

This argument is specific to full14-pair frames. It does not classify all lattices, assume integrality of Λ, or prove global rigidity of arbitrary near-extremal sections. The finite-group and sign-averaging steps are mathematical arguments stated here; the replay checks their finite arithmetic inputs, not a formal Lean proof of the entire theorem.
