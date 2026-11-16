# 15D near-frame arithmetic gap

Assume minimum length 1 and original density bound Δ≤0.0256118267. If 15 opposite vector pairs share one nonzero mod-2 class, write S=Σ(|v_i|²−2). Then

**ρ>1024√2 implies S>1/5.**

For contradiction assume S≤1/5. The signed near-frame determinant lemma (n≥9, S≤8/(n−1)) gives det T≥1 for T=[v_i]/√2. The midpoint constraints give metric distortion D≤1+(1+√15)S/4≤199/160. The sublattice T(D15/√2) has index I=ρ 2^(−13/2)det T. Hadamard and π>3.14159 give exactly

I² ≤ 2(0.0256118267)²(2027025)²/(3.14159)^14 ×(151/150)^15 < 676.

Thus I≤25. Exact canonical histogram enumeration gives the following maximum, over order-p points of R15/D15, of the least physical squared norm among nonzero multiples:

p: 5, 7, 11, 13, 17, 19, 23
bound: 4/5, 36/49, 88/121, 108/169, 178/289, 228/361, 323/529.

Each is at most4/5 and D·4/5≤199/200<1, excluding these primes from I. The portable verifier verify_arithmetic.py independently enumerates all histograms for primes5,7,11,13 and uses exact parity-support average bounds for17,19,23. Its raw distance is Σr_i²/p² plus, if centered integer residues sum odd, 1−2max|r_i|/p. Physical squared norm is half this.

Additionally the portable verify_arithmetic.py enumerates all26720 exact-order6 histograms and gives physical witness≤7/9. This also contradicts D≤199/160. Completeness: any a/6 with a∈D15 can be reduced coordinatewise by shifts in2Z and sign changes to numerators in{0,...,6}; these preserve D15. The script tracks the parity of nearest-integer reductions, including half-coordinate ties, and discards a histogram only if a proper nonzero multiple lies in D15. Thus parity-kernel lifts are included.

Finite abelian group Λ/T(D15/√2) cannot have order divisible by6, since such a group contains an element of order6. Its order is≤25 and has no prime≥5. The maximum remaining order is16. Therefore ρ≤16·2^(13/2)/det T≤1024√2.

All numerical comparisons in the proof are rational. The source density bound is an input; this file does not re-certify the published bound.

## Attribution and portable dependencies

The determinant step is an application of Géraud-Stewart–Naccache, *A Focus on Hermite Constant11: the certificates and macadamia0052*, AppendixD, Lemma32, pp44–45, https://eprint.iacr.asia/2026/2826.pdf . The accompanying determinant proof/checker is an independent reconstruction, not a novelty claim for that lemma.

Run `python run_arithmetic.py` from any directory. It replays the bundled original Zubrilina15D prior using `../prior/taylor.py` and `prior_certificate.json`, verifies the determinant identities, and recomputes all arithmetic/torsion cases. The upstream prior scripts/data are pinned by the package provenance/manifest; do not substitute floating-point optimization logs for these inputs. No network, C++ compiler, absolute workspace path, or cached PASS report is needed for this replay. Python3 and SymPy are required.

The two enumerations differ intentionally: the original full histogram screen gave sharper witnesses for primes17,19,23; the portable independent verifier uses coarser exact parity-support averages, which still satisfy the required strict inequality after multiplying by199/160. Both prove the same theorem.
