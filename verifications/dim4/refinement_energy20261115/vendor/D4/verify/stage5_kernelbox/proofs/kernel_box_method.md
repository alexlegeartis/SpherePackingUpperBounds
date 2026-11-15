# Kernel-box certificate: -K >= lambda_j q(t) on B_1 and B_5

Family: all packings in R^4, through the raw V7 ledger and the reallocation of `conditional_theorem.md`. Evidence: CERTIFIED by interval arithmetic (python-flint, Arb balls at 256 bits) with a written soundness argument below.

## 1. Claim

Let a = 133/250 and q(t) = t^2 (t^2 - a^2). Let K be the raw V7 pair kernel, the polarization of `src/raw_v7.py` (as in `src/kernel_poly.py`). For j = 1, 5 put

    B_j = {(y, z) : |y|, |z| in [1, R_j], t = cos(y, z) in [-1, -a]},

with R_1 = 501473/500000 and R_5 = 102597/100000.

- -K(y, z) >= (3/1000000) q(t) on B_1.
- -K(y, z) >= (1/1000000) q(t) on B_5.

K depends on y, z only through r = |y|, s = |z| and t. So each claim is a statement on the box [1, R_j]^2 x [-1, -a] in (r, s, t).

## 2. An explicit formula for K

Because the anchored energy is a quadratic form in the star, K is the sum of the cross terms of each component of the raw V7 energy. `src/kformula.py` writes these out, with Z_d = U_d/(d + 1):

- Gram part: sum over d = 0..18 of Z_d(t) times sum over rays of (w/DEN) alpha(r) alpha(s), with the Bernstein features alpha of `data/v7/kernel.json`. Degree 0 is the origin block E0.
- Spherical rays: -weight g(t) when both radii are at most the row radius.
- Count and conditional polynomials: their polarizations. These are constants on each radius class.
- Gauge: (G(s) q0(r)/(1 - psi_r) + G(r) q0(s)/(1 - psi_s))/48, where q0 is the unlifted pair of raw_v7 and G, psi are the gauge indicator and weight.
- Pos, the part that enters -K with a plus sign: caps (weight/2)(rw(s) P(y|z) + rw(r) P(z|y)), tents (weight/2)(rw(s) T(y|z) + rw(r) T(z|y)), allocations (weight/2)(rw(r) + rw(s)) cap(|y + z|^2), and graphs (the weight per qualifying pair). All weights and radial weights rw are nonnegative, and P, T and cap lie in [0, 1].

So -K = smooth + Pos, where smooth = -(Kquad + Ksph + Kcount + Kcond + Kgauge).

Independent check. `src/validate_kformula.py` compares kformula with K computed by polarization of the raw_v7 evaluator (`src/kernel_poly.py`, which calls `src/raw_v7.py` on 1- and 2-point stars with exact rationals). It uses 72 rational points over nine radius pairs, including every break R_1, ..., R_5 and 103/100, with t from -0.997 to -0.23. It also compares q0 with `raw_v7.unlifted_pair` at 9 radii. Every kformula ball must meet the raw_v7 enclosure. The derivation, not the point check, identifies the two formulas. The point check guards against transcription errors. raw_v7 is written separately from kformula: it works with exact rationals and outward rational square-root bounds, and it evaluates whole stars rather than the cross-term formula.

A second check, `src/coord_check_kernelbox.py`, uses raw_v7 only. It checks that E on random rational 3- and 4-point stars equals A + sum of singleton terms + sum of 2K over pairs (additivity), and that -K - (3/1000000) q >= 0 at 88 rational points of B_1 near the dip and at band corners.

## 3. Box enclosures and kinks

`src/kcert.py` bounds -K - lambda q from below on a box by a second-order centered form. Each quantity carries:

- its value and gradient at the box centre (Arb balls), and
- enclosures of its value, gradient and Hessian over the whole box.

On the box, f lies in c + g . d + (1/2) d^T H d with |d_i| <= half-width_i (Taylor's theorem with Lagrange remainder, H enclosing the Hessian over the box). This is intersected with the plain value enclosure. Every operation used (sums, products, 1/x away from 0, exp, sqrt away from 0) is C^2 on the box, and the code asserts the side conditions.

The kinks are the cap indicators and clipped probabilities, the tent relu, the allocation step capacities and the graph indicator. On a box where a kink is decided, kcert uses the exact branch. Where it is undecided, kcert replaces that term by a constant lower bound of the term over the box (zero derivatives). For example, relu(x) >= 0, clip(x, 0, 1) >= max(0, inf x), an indicator [L <= 0] >= 0, and a nonincreasing step capacity is at least its value at sup |y + z|^2.

This is sound because -K is nondecreasing in each kinked term:

- In Pos each kinked term has a nonnegative coefficient.
- In q0 the self terms enter as +weight rw (M - P), +weight rw (capacity - relu) and +weight rw (1 - cap), and q0 enters -K with coefficient -(1/48)/(1 - psi) <= 0. So -K is again nondecreasing in P, relu and cap.
- On B_5 the self terms are all decided. kcert asserts this ("self cap undecided" would abort).

So on every box kcert computes the centered form of a C^2 function that is at most -K - lambda q. A nonnegative lower bound proves the claim on that box.

Radius breaks. The features are piecewise on the radial regions [R_k, R_(k+1)] with a closed-inner convention: the face r = R_k belongs to region k - 1. A root box whose lower radius is a break R_k > 1 is marked open on that side. Threshold decisions on it use lo + 10^-40, below every gap between distinct radius thresholds. The face itself is covered by the neighbouring lower box.

## 4. The subdivision and its replay

`src/kbox_run.py` tiles B_j by root boxes: each pair of radius regions [R_k, R_(k+1)] (k < j), times 94 equal t-chunks of [-1, -a]. It bisects each root at midpoints until every leaf has a nonnegative lower bound, and the run fails if a box still has a negative bound when none of its sides is wider than 2 * 10^-7. It saves the split tree of every root.

| j | lambda_j | roots | leaves | smallest leaf lower bound |
|---|---|---|---|---|
| 1 | 3/10^6 | 94 | 286 | 2.2e-8 |
| 5 | 1/10^6 | 2350 | 4663 | 1.6e-9 |

The trees are `data/kernel-box/kbox_j1.json` and `kbox_j5.json`. They store only the split structure, never leaf values.

`src/check_kernelbox.py` gate 4 replays both trees without search. The root list must equal the tiling of B_j computed by `kbox_run.roots`, every split must be a midpoint bisection on dimension 0, 1 or 2, and kcert must give a nonnegative lower bound at every leaf, recomputed from scratch.

The binding point is the dip at radii (1, 1) and t about -0.864. There -K = 1.0773e-6 and -K/q = 3.114e-6, so lambda_1 = 3/10^6 is 96% of the pointwise optimum. lambda_1 = 4/10^6 is false there, and both the replay and the search reject it (negative control).

## 5. From the kernel bound to the ledger

`far_pair_lemma.md` turns the two kernel bounds into H_R >= lambda_j max(Q_N, 0) for every census group (N, j). `src/ledger.py` then computes

    a_new = min over the 124 groups of floor_g + lambda(j) max(Q_N, 0),

with lambda(1) = 3/10^6 and lambda(j) = 1/10^6 for j = 2..5. The minimum is attained only at (24, 1): a_new = a_plus + (3/10^6) Q_24 = a_plus + 65799/1562500000 = 0.0038258783318289675... The next group is (26, 5) at 0.0038387440017070..., a gap of 1.29e-5.

    U_new = U_raw - a_new/8 - 3/8000000 = 0.629521681652091506... < 0.629521681653.

The previous ceiling was 0.629526945573, so the improvement is 65799/12500000000, about 5.26e-6.

## 6. Hypotheses inherited from other stages

1. Stage 1: the base certificate signs, including the cyclic sign Gamma <= 0, and the definition of K by raw_v7.
2. Stages 2 and 3: the 124 census floors and the group definition (j is the first near radius R_j that contains all N near neighbours).
3. Stage 4: the graph saving 3/8000000.
4. The density theorem of `conditional_theorem.md`.
