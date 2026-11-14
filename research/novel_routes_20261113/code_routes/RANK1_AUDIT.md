# Independent audit of the 71D rank-one construction

The higher-parent search found a primitive vector v of squared norm 10 in
Nebe's extremal 72D even unimodular lattice. Taking the two adjacent labels
0 and1 in the projected dual line gives a valid 71D antipode packing with

    delta = 2/sqrt(10) * (79/40)^(71/2)
          = 19663555319.937... .

This note independently checks that this is the largest possible value of
the rank-one antipode formula based solely on the parent minimum 8. This
does not establish a world record among all 71D packings.

## Exhausting the one-dimensional label choices analytically

For a primitive normal of squared norm a, its projected dual line has
spacing 1/sqrt(a). Any set of s distinct labels has squared diameter at
least (s-1)^2/a; equality holds for consecutive labels. Its antipode formula
is therefore at most

    D_s(a) = s/sqrt(a) * ((8-(s-1)^2/a)/4)^(71/2),

whenever the bracket is positive. Since the parent is even and has
minimum 8, a is an even integer at least 8.

For s=1, the largest value is the bare norm 8 section, namely 2^34.

For s=2, differentiating the logarithm of D_2(a)^2 gives a quantity with
the sign of72-8a. Thus the continuous maximum is at a=9, and the even
integer maximum must be at 8 or 10. The exact squared ratio is

    D_2(10)^2 / D_2(8)^2 = (4/5)*(316/315)^71 > 1.

For s>=3, put b=a/(s-1)^2. Then

    D_s(a) = [s/(s-1)] / sqrt(b) * ((8-1/b)/4)^(71/2).

The continuous maximum in b is again at 9. Also s/(s-1)<=3/2, so the
whole s>=3 family is bounded by

    (1/2)*(71/36)^(71/2).

The exact squared ratio of the proposed packing to that upper bound is

    (8/5)*(711/710)^71 > 1.

The norm 10 construction also beats the s=1 bound. These comparisons
prove the claimed formula optimum, without relying on a finite scan of
normal norms or cluster sizes. Dropping possible restrictions from label
distinctness modulo the normal lattice only enlarges the relaxation used
for the upper bound, so it does not undermine the proof.

## Arithmetic audit of the higher-parent verifier

The norm 10 normal has gcd 1 in the published ambient basis, so it is
primitive. Its Gram product is exactly 10, and the two labels are 0,1.
The source verifier additionally checks that the ambient Gram agrees with
the archived catalogue entry, is even, positive definite and determinant 1.
The minimum 8 theorem is inherited from Nebe's published construction.

The higher-parent bare-section comparison constants in ranks 1 through 8
also check out: they follow from the exact Hermite constants in those
ranks, applied to the normal lattice of minimum at least 8. They bound
the density *formula that retains the parent's minimum 8*. They are not
upper bounds for all sections if a section's actual minimum is larger.

Replay:

    python research/novel_routes_20261113/code_routes/verify_rank1_antipode.py

## Priority and benchmark status

The retrieved current Cohn table covers dimensions 1–48, 56, 64, 72. It has
no 71D entry. The Nebe–Sloane table is likewise sparse above 48. Targeted
searches found the norm 8 orthogonal section of the putative/extremal 72D
lattice recorded in OEIS A064163 since 2001, but did not locate a current
comprehensive density table in dimensions 65–71 or a prior printed value
for the norm 10 antipode packing. Absence of such a source is not evidence
that no denser construction is known.

The candidate is 14.456956...% denser than the bare norm 8 section and
about 0.09486% denser than the two-label norm 8 antipode formula. These are
checked comparisons with stated constructions, not a substitute for a
world-record comparison. Priority remains unestablished.

Sources:

- Cohn's current sphere-packing table: <https://cohn.mit.edu/sphere-packing/>.
- Conway–Sloane, *The Antipode Construction for Sphere Packings*:
  <https://neilsloane.com/doc/anti.pdf>.
- The recorded 71D norm 8 section:
  <https://oeis.org/A064163>.
- The published parent and its relatives:
  <https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/ext72.html>.
