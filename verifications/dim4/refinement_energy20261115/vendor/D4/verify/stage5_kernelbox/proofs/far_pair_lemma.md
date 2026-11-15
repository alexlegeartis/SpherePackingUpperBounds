# Far-pair lemma: H_R >= lambda_j max(Q_N, 0)

PROVED by written argument, given the kernel-box certificates (H3). Scope: every anchor of every packing in R^4, in the raw V7 local-energy ledger with the reallocation K_R = K + H_R (`conditional_theorem.md`). This is a lower bound on H_R, not a density bound by itself.

## Setting

- K is the pair kernel of the base certificate, the polarization of `src/raw_v7.py`. K(y, z) = 0 when |y| or |z| exceeds 7/4.
- For points y, z of the 7/4 star of an anchor x (placed at the origin), Gamma(y, z) = K(y, z) + K(-y, z - y) + K(-z, y - z).
- H_R,x = sum over ordered pairs y != z of the star of -Gamma(y, z)/k_R(y, z), where k_R is the number of triangle vertices whose two incident sides are at most R = 7/4. Terms with Gamma = 0 are omitted, and k_R >= 1 whenever Gamma != 0.
- R_1, ..., R_5 are the kernel breaks 501473/500000, 1005643/1000000, 50379/50000, 50767/50000 and 102597/100000 (`data/v7/kernel.json`). The cap of R_j is s_j = 1 - 1/(2 R_j^2).
- a = 133/250 and q(t) = t^2 (t^2 - a^2).
- B_j is the set of pairs (y, z) with |y|, |z| in [1, R_j] and cos(y, z) in [-1, -a].

## Hypotheses

1. (H1) Cyclic sign. Gamma <= 0 on every admissible triangle. This is the accepted cyclic sign of the base certificate (stage 1), kept by the reallocation.
2. (H2) Allocation rule. H_R is the reallocation of `conditional_theorem.md`.
3. (H3) Kernel bound. -K(y, z) >= lambda_j q(t) on B_j, with lambda_1 = 3/1000000 on B_1 and lambda_5 = 1/1000000 on B_5. These are the kernel-box certificates (`kernel_box_method.md`), replayed leaf by leaf in this stage. Since B_j lies in B_5 for j <= 5, lambda_5 serves j = 2, ..., 5.

## Statement

Let x be an anchor and j in {1, ..., 5}. Let y_1, ..., y_N be N distinct packing points with 1 <= |y_i - x| <= R_j, placed relative to x. Put lambda(1) = lambda_1 and lambda(j) = lambda_5 for j >= 2. Then

    H_R,x >= lambda(j) max(Q_N, 0),   Q_N = N^2 (1/8 - a^2/4) - N (1 - a^2) = N (13561 N - 179244)/250000.

Q_N > 0 exactly when N >= 14.

| N | Q_N (exact) | Q_N (decimal) |
|---|---|---|
| 24 | 43866/3125 | 14.037120 |
| 25 | 159781/10000 | 15.978100 |
| 26 | 1126723/62500 | 18.027568 |

## Proof

Write u_i = y_i/|y_i| and t_ik = u_i . u_k.

Step 1. Positive type. On S^3 the Gegenbauer parameter is (4 - 2)/2 = 1, and C_k^1 = U_k. With U_4 = 16t^4 - 12t^2 + 1, U_2 = 4t^2 - 1 and U_0 = 1,

    q = (1/16) U_4 + (3/16 - a^2/4) U_2 + (1/8 - a^2/4) U_0
      = (1/16) U_4 + (14593/125000) U_2 + (13561/250000) U_0.

All three coefficients are positive. Each U_k(u . v) is a positive definite kernel on S^3 (Schoenberg 1942). So the sum over all ordered pairs (i, k), diagonal included, of U_k(t_ik) is >= 0 for k >= 1, and it equals N^2 for k = 0. The full sum of q(t_ik) is therefore at least (13561/250000) N^2. The diagonal contributes N q(1) = N (1 - a^2) = 44811 N/62500. So

    sum over i != k of q(t_ik) >= Q_N.        (1)

Step 2. Cap. Take i != k, r = |y_i| and rho = |y_k|, both in [1, R_j]. Unit separation gives t_ik <= h(r, rho) = (r^2 + rho^2 - 1)/(2 r rho). For fixed rho >= 1, h = r/(2 rho) + (rho^2 - 1)/(2 rho r) is convex in r > 0, and likewise in rho. So the maximum over [1, R_j]^2 is at a corner. The corner values are h(1, 1) = 1/2, h(1, R) = R/2 and h(R, R) = s_j. Also R/2 <= s_j is equivalent to (R - 1)(R^2 - R - 1) <= 0, which holds for 1 <= R <= (1 + sqrt 5)/2, and every R_j is in that range. Hence t_ik <= s_j. Exactly, s_5 = 0.5249922... < a = 0.532, and s_j <= s_5.

Step 3. Sign split. q(t) = t^2 (t^2 - a^2) <= 0 for |t| <= a. By Step 2 every pair has t_ik <= s_j <= a, so every pair with t_ik >= -a has q(t_ik) <= 0. Let F be the set of ordered pairs i != k with t_ik < -a. Then (1) gives

    sum over (i, k) in F of q(t_ik) >= sum over i != k of q(t_ik) >= Q_N.        (2)

Step 4. Pairs in F are far. For (i, k) in F with r, rho >= 1,

    |y_i - y_k|^2 = r^2 + rho^2 - 2 r rho t_ik > r^2 + rho^2 + 2a r rho >= 2 + 2a = 383/125 > 49/16.

So |y_i - y_k| > 7/4. In the triangle (x, y_i, y_k) only the vertex x has both incident sides at most 7/4, so k_R = 1. The vectors y_k - y_i and y_i - y_k are longer than 7/4, so K(-y_i, y_k - y_i) = K(-y_k, y_i - y_k) = 0 and Gamma(y_i, y_k) = K(y_i, y_k).

Step 5. Assembly.

- By (H1) every term -Gamma/k_R of H_R,x is >= 0, so dropping terms only lowers the sum.
- The points y_i lie in the 7/4 star, since R_j < 7/4. Keep only the ordered pairs in F.
- By Step 4, each kept term equals -K(y_i, y_k).
- Each (y_i, y_k) with (i, k) in F lies in B_j, so (H3) gives -K >= lambda(j) q(t_ik).
- With lambda(j) >= 0 and (2),

      H_R,x >= sum over F of -K(y_i, y_k) >= lambda(j) sum over F of q(t_ik) >= lambda(j) Q_N.

Also H_R,x >= 0 by (H1). Hence H_R,x >= lambda(j) max(Q_N, 0). QED.

## Remarks

- No localization. The bound does not need the near code to be close to the 24-cell or to be antipodal. It uses only the count N, the cap and the sign of q.
- Group use. For a census group (N, j), the group's N near neighbours are the y_i. They lie in [1, R_j] by the group definition (first containing radius R_j).
- Dropped terms. Pairs with t >= -a and all non-near pairs of the star are dropped. Their terms are >= 0 by (H1).
- Only the certified lambdas are used. The pointwise optimum of -K/q on B_1 is about 3.114e-6, at radii (1, 1) and cosine about -0.864.
- Range of j. The statement needs R_j <= 1.618 for Step 2 and 2 + 2a > 49/16 for Step 4. It is stated only for j <= 5.

## Checks

`src/hrlemma_check.py` checks with exact rationals: the Gegenbauer identity and the coefficient signs, q(1), the cap and corner inequalities for j = 1..5, 383/125 > 49/16, the Q_N formula, the threshold N >= 14 and Q_24, Q_25, Q_26. Its last check only reads the summary block of the two trees. The trees themselves are replayed leaf by leaf in `src/check_kernelbox.py` gate 4. The sign of q on [-a, a] comes from the factorization in Step 3.
