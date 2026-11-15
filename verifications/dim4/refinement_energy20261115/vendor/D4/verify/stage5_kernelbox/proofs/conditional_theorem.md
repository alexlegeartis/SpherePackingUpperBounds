# The reallocated kernel K_R = K + H_R and the density theorem

Family: all packings in R^4. Source: the D4 full-support allocation of the project, cleaned for this package. Proof labels: PROVED (written argument) unless marked otherwise.

## Notation

- f(x) = exp(-2|x|^2) F(2|x|^2) is the scalar function of the base certificate (`data/v7/scalar.json`). Its Fourier transform is nonnegative, and I = fhat(0) = pi^2/4.
- R = 7/4 is the support radius.
- K is the pair kernel of the base certificate. It is the polarization of the anchored energy E evaluated by `src/raw_v7.py`: for a != b, K(a, b) = (E({a, b}) - E({a}) - E({b}) + A)/2, where A = E(empty) = 8 U_raw - f(0). K(a, b) = 0 if |a| or |b| exceeds 7/4.
- For a triangle (0, y, z), the cyclic term is Gamma(y, z) = K(y, z) + K(-y, z - y) + K(-z, y - z). It has the same value at all three anchors of the triangle.
- k_R(y, z) is the number of triangle vertices whose two incident sides both have length at most R. It is 0, 1 or 3.
- C(U) = 8U - f(0).

## 1. The density theorem

Hypothesis H(U). For every finite unit-separated set P in R^4,

    S(P) = sum over ordered pairs x != y in P with |x - y| <= R of f(x - y) <= C(U) |P|.

Theorem. H(U) implies that every packing of unit balls (centres at mutual distance at least 1) in R^4 has upper density at most U. The density here is (pi^2/32) times the centre density. For a proved constant C the ceiling is U = (f(0) + C)/8.

Proof.

(a) Sign beyond R. For |x| > R, K(x, 0) = K(x, x) = 0 because K is supported on radii in [1, R]. The accepted pair sign f + 2K(x, 0) + K(x, x) <= 0 then gives f(x) <= 0.

(b) Finite inequality. For a finite unit-separated P with N points,

    sum over x, y in P of f(x - y) = N f(0) + S(P) + sum over |x - y| > R of f(x - y) <= N (f(0) + C).

(c) Periodic packings. Let Lambda + {p_1, ..., p_N} be periodic with covolume V. Apply (b) to a patch of M^4 periods. Points farther than R from the patch boundary have all their R-neighbours inside. The positive part of f has support |d| <= R and f is bounded, so the boundary share is O(M^3). Per period, sum over i and over packing points y of f(p_i - y) <= N (f(0) + C); the terms with |d| > R are at most 0 and only help. Poisson summation with fhat >= 0 bounds the left side below by N^2 I / V. So N/V <= (f(0) + C)/I, and the density is (pi^2/32) N/V <= (pi^2/32)(f(0) + C)/(pi^2/4) = (f(0) + C)/8.

(d) Arbitrary packings. Periodic packings approximate the upper density of any packing. The finite-patch argument of the local-energy theorem starts from (b) with A replaced by C and gives the same bound.

## 2. Where the reallocation enters, with no double count

Section 1 uses only f. The reallocated kernel K_R = K + H_R explains what H(U) means and proves the instance used in stage 5.

- H_R(y, z) = -Gamma(y, z)/k_R(y, z), summed at an anchor over ordered pairs y != z of its star with |y|, |z| <= 7/4. Terms with Gamma = 0 are omitted. Gamma is nonzero only if some vertex has both incident sides in [1, R], so k_R >= 1 whenever Gamma != 0.
- H_R vanishes on the diagonal and at the origin. So A = K_R(0, 0) = K(0, 0), and the pair term f + 2K_R(x, 0) + K_R(x, x) equals the accepted pair term. The accepted pair sign is unchanged.
- The k_R anchors of a triangle each receive -Gamma/k_R, so the three allocations sum to -Gamma. The new cyclic sum is 0 <= 0, so the accepted cyclic sign holds.
- The accepted cyclic sign Gamma <= 0 gives H_R >= 0. So E_R = E + H_R >= E for every admissible neighbourhood.
- K_R is symmetric, bounded, Borel and supported in [1, R].
- Grouping the energy by anchors, with zero cyclic sum and the completed pair identically zero on [1, R], gives sum over x of E_R,x = A N - S(P).
- So H(U) is equivalent to the averaged floor (1/N) sum over x of E_R,x >= A - C(U) = 8 (U_raw - U).

Consequence. If every anchor has E_R,x >= a, then H(U) holds with C = A - a, and the ceiling is U = U_raw - a/8.

No double count. The ceiling formula contains only f(0) and C. The floor a appears once. The graph charges (4769/200000000 in U_raw) sit inside A and K and never appear separately. The broad graph saving of stage 4 raises the census floors by 3/1000000 and so enters a once.

## 3. How stage 5 uses it

Stages 2 and 3 give, for each of the 124 census groups g = (N, j), a floor floor_g for E at an anchor of group g (`census_floors.json`, field `new_floor`). Stage 4 raises it by the graph slack 3/1000000. The far-pair lemma (`far_pair_lemma.md`) gives H_R >= lambda_j max(Q_N, 0) at the same anchor. So

    E_R,x = E_x + H_R,x >= floor_g + lambda_j max(Q_N, 0) + 3/1000000,

and the minimum over the 124 groups is a_new + 3/1000000. The ceiling is

    U_new = U_raw - a_new/8 - 3/8000000 < 0.629521681653.
