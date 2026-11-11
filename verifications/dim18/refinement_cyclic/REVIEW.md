# Independent cyclic certificate review

The frozen witness is dimension 18, with five radial hats in squared radius, solid-harmonic angular degrees 0 through 16, and Gaussian polynomial degree 31. Its exact exclusion and truncation squares are 9/4 and 26999999999999997/10000000000000000. Certificate SHA256: 571075177f161a68a671cc6b9057cfb6d7f8f8c7613853698c4a164844fbfa4e.

## Mathematical implication

The printed Theorem 1.1 in Cohn–de Laat–Salmon (https://arxiv.org/html/2206.15373v1) assumes individual kernel nonpositivity on the admissible annular pair domain. This witness instead uses a cyclic generalization. It must not be described as satisfying that individual-sign hypothesis without a separate check. The following finite-sum argument supplies the required implication.

Let f be the unnormalized Gaussian polynomial and K the unnormalized supported PSD kernel. Normalize packing separation to r=sqrt(9/4). At each packing center p, sum K over the displacement set consisting of zero and its neighbors within R=sqrt(26999999999999997/10000000000000000). PSD gives a nonnegative sum. Its expansion is K(0,0), the sum over neighbors of 2K(x,0)+K(x,x), and the ordered distinct-neighbor terms.

Sum this expansion over the N center representatives of a periodic packing. Enlarge the period cell so no short displacement or short triangle is identified incorrectly. Ordered distinct-neighbor terms regroup as twice the cyclic sums over translation classes of unordered triangles. Each such cyclic sum is nonpositive: if all three sides are at most R, use the compact cyclic certificate; with one side greater than R only its opposite vertex term survives and the far certificate applies; with at least two long sides every term vanishes. Consequently the summed pair expression is at least -N K(0,0). The verified pair inequality f(x)+2K(x,0)+K(x,x)<=0 and physical tail f<=0 imply the off-diagonal sum of f is at most N K(0,0).

Because f is Schwartz and its Fourier transform is nonnegative with H0=fhat(0)>0, periodic Poisson summation yields N^2 H0/V <= N(f(0)+K(0,0)). Thus density is at most vol(B_(r/2)^18)*(f(0)+K(0,0))/H0. Periodic approximation of finite patches, with separating boundary gaps whose volume fraction tends to zero, gives the unrestricted packing bound. No lattice-only, average-neighbor-count, or integer-count assumption is used.

The solid-harmonic factor is (sqrt(u*v)/r^2)^ell P_ell^(18)((u+v-w)/(2sqrt(u*v))). Its positive radial factors can be absorbed into spherical harmonic features, preserving PSD. Hats and all angular terms vanish outside the closed annulus; the origin has a separate degree-zero coordinate. Continuity of this kernel is not needed for the finite sums.

## Domain geometry

Cyclic cell indices i<=j<=k suffice because the cyclic polynomial is symmetric in the three squared side lengths. The full Cartesian product is checked even where it contains unordered or unrealizable triples; this is stronger than needed. At repeated cell indices there is no ordering restriction within the cell. Every actual triangle is included after permuting side labels.

For exactly one long side, its squared length w satisfies R^2<w<=(sqrt(u)+sqrt(v))^2<=4R^2. Checking K on the complete rectangular range w in [R^2,4R^2] is therefore sufficient. Our independent certificate checks that entire range, including unrealizable triangles, with no geometric exclusions. Closed boundaries are included in both checks.

## Independent computation

algebra_pair.py independently checks 17 exact PSD blocks (86 positive rational pivots), the direct monomial Fourier transform, global Fourier positivity and physical tail via rational Sturm sequences, and the complete solid-harmonic pair inequality using independently reconstructed Arb polynomial Taylor enclosures. The pair check passes in 16 leaves of depth at most 6. A rational Machin upper bound for pi proves the final squared density is less than 0.013284728 squared. Arb at 384 bits gives density [0.01328472704022940068952086716325760124601533840031853772476992896674363378012149103489871385717865691637676166429868 +/- 7.53e-117].

domains.py builds solid harmonics from their explicit rising-factorial formula rather than the producer's recurrence. It separately performs exact affine monomial substitution and exact power-to-Bernstein conversion. Every subsequent midpoint de Casteljau subdivision is rational. Floating point only chooses which axis to split; it never proves a sign or excludes a region. All 20 compact cyclic boxes and all 10 far boxes must pass before the script prints PASS. See domains.log and domains_audit.json for the final run status.

The frozen audit source hash above identifies the original export. Registration reconciles only metadata: strict upper 0.013284728, completed audit status and solid-harmonic family labels. All rational polynomials, matrices, radii, degrees, normalization and objective are unchanged by exact JSON comparison. Registered certificate SHA256: 84fe1f8779e84c4bb7f4643ad7afc9afcc9761dd9f793fd0f33e934319055b82. Both paths are rerun against this file.

## Benchmark and limitations

A fresh search on 11 November 2026 retrieved Cohn's indexed general-packing table https://cohn.mit.edu/sphere-packing/, which lists upper 0.0134853404450862 in dimension 18 and attributes it to Afkhami-Jeddi et al., https://arxiv.org/abs/2006.02560. The newer 2512.04829 paper covers dimensions 4–7 and 9–16, not 18. Refreshed searches found no stronger finite 18D general bound. The available Cohn index is four months old, so this is a pinned-source comparison, not a worldwide-priority certification.

This is a computer-assisted proof audit, not proof-assistant formalization. It trusts rational arithmetic and Arb, standard Fourier identities and the Gegenbauer addition theorem. It makes no novelty claim for the cyclic method.

The prior newly certified repository bound in dimension18 was 0.013408931; the new cyclic witness also improves that result. The exact truncation square is the recorded rational 2.6999999999999997, not rounded to 27/10.
