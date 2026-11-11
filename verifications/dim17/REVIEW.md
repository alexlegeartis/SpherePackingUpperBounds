# Independent cyclic certificate review

The frozen witness is dimension 17, with five radial hats in squared radius, solid-harmonic angular degrees 0 through 16, and Gaussian polynomial degree 31. Its exact exclusion and truncation squares are 17/8 and 51/20. Certificate SHA256: 6313cfe172ca945cef7e16f4808bf2457229575dc9883a631a0e8b406ad9c3ea.

## Mathematical implication

The printed Theorem 1.1 in Cohn–de Laat–Salmon (https://arxiv.org/html/2206.15373v1) assumes individual kernel nonpositivity on the admissible annular pair domain. This witness instead uses a cyclic generalization. It must not be described as satisfying that individual-sign hypothesis without a separate check. The following finite-sum argument supplies the required implication.

Let f be the unnormalized Gaussian polynomial and K the unnormalized supported PSD kernel. Normalize packing separation to r=sqrt(17/8). At each packing center p, sum K over the displacement set consisting of zero and its neighbors within R=sqrt(51/20). PSD gives a nonnegative sum. Its expansion is K(0,0), the sum over neighbors of 2K(x,0)+K(x,x), and the ordered distinct-neighbor terms.

Sum this expansion over the N center representatives of a periodic packing. Enlarge the period cell so no short displacement or short triangle is identified incorrectly. Ordered distinct-neighbor terms regroup as twice the cyclic sums over translation classes of unordered triangles. Each such cyclic sum is nonpositive: if all three sides are at most R, use the compact cyclic certificate; with one side greater than R only its opposite vertex term survives and the far certificate applies; with at least two long sides every term vanishes. Consequently the summed pair expression is at least -N K(0,0). The verified pair inequality f(x)+2K(x,0)+K(x,x)<=0 and physical tail f<=0 imply the off-diagonal sum of f is at most N K(0,0).

Because f is Schwartz and its Fourier transform is nonnegative with H0=fhat(0)>0, periodic Poisson summation yields N^2 H0/V <= N(f(0)+K(0,0)). Thus density is at most vol(B_(r/2)^17)*(f(0)+K(0,0))/H0. Periodic approximation of finite patches, with separating boundary gaps whose volume fraction tends to zero, gives the unrestricted packing bound. No lattice-only, average-neighbor-count, or integer-count assumption is used.

The solid-harmonic factor is (sqrt(u*v)/r^2)^ell P_ell^(17)((u+v-w)/(2sqrt(u*v))). Its positive radial factors can be absorbed into spherical harmonic features, preserving PSD. Hats and all angular terms vanish outside the closed annulus; the origin has a separate degree-zero coordinate. Continuity of this kernel is not needed for the finite sums.

## Domain geometry

Cyclic cell indices i<=j<=k suffice because the cyclic polynomial is symmetric in the three squared side lengths. The full Cartesian product is checked even where it contains unordered or unrealizable triples; this is stronger than needed. At repeated cell indices there is no ordering restriction within the cell. Every actual triangle is included after permuting side labels.

For exactly one long side, its squared length w satisfies R^2<w<=(sqrt(u)+sqrt(v))^2<=4R^2. Checking K on the complete rectangular range w in [R^2,4R^2] is therefore sufficient. Our independent certificate checks that entire range, including unrealizable triangles, with no geometric exclusions. Closed boundaries are included in both checks.

## Independent computation

algebra_pair.py independently checks 17 exact PSD blocks (86 positive rational pivots), the direct monomial Fourier transform, global Fourier positivity and physical tail via rational Sturm sequences, and the complete solid-harmonic pair inequality using independently reconstructed Arb polynomial Taylor enclosures. The pair check passes in 14 leaves of depth at most 5. A rational Machin upper bound for pi proves the final squared density is less than 0.018335552 squared. Arb at 384 bits gives density 0.01833555157029049139332536094236693008327772532524438023943130167025927979684206224099336131303752154313179742348101, with enclosure radius below 7.5e-117.

domains.py builds solid harmonics from their explicit rising-factorial formula rather than the producer's recurrence. It separately performs exact affine monomial substitution and exact power-to-Bernstein conversion. Every subsequent midpoint de Casteljau subdivision is rational. Floating point only chooses which axis to split; it never proves a sign or excludes a region. All 20 compact cyclic boxes and all 10 far boxes must pass before the script prints PASS. See domains.log and domains_audit.json for the final run status.

The frozen audit source SHA256 above identifies the original discovery export. Registration changes only metadata: the strict upper is 0.018335552, the status records the completed audits, and the kernel-family labels identify the solid-harmonic cyclic representation. All rational polynomials, matrices, radii, degrees, normalization and objective are unchanged by an exact JSON comparison. Registered certificate SHA256: 32a9c22f599c91d80d35a585cbe8ca2e40d53a189d23d6f7ba8c400835f307e4. Both verification paths are rerun against this registered file.

## Benchmark and limitations

A fresh search on 11 November 2026 retrieved Cohn's indexed general-packing table https://cohn.mit.edu/sphere-packing/, which lists upper 0.0184640903350649 in dimension 17 and attributes it to Afkhami-Jeddi et al., https://arxiv.org/abs/2006.02560. The newer 2512.04829 paper covers dimensions 4–7 and 9–16, not 17. Refreshed searches found no stronger finite 17D general bound. The available Cohn index is four months old, so this is a pinned-source comparison, not a worldwide-priority certification.

This is a computer-assisted proof audit, not proof-assistant formalization. It trusts rational arithmetic and Arb, standard Fourier identities and the Gegenbauer addition theorem. It makes no novelty claim for the cyclic method.
