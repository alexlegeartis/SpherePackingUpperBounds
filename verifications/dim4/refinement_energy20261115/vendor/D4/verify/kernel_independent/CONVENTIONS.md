# Conventions adopted in kernel.py

`kernel.py` is written from `data/v7`, Section 3 of `paper.tex`,
`verify/stage1_v7/COMPONENTS.md` and the proof notes in
`verify/stage1_v7/proofs/`. Where those leave a choice open, I read the stage 1
code to learn its convention and adopted it. Each such choice is listed here
with the stage 1 file and line. Paths below are relative to
`verify/stage1_v7/record/current/`. "Coordinator" means
`kernel_integration_route/verify_combined_binary_six_exterior_graphs_small_cap_bounds.py`.
No stage 1 code is imported by `kernel.py`.

## Domain, pieces and break points

1. **Pieces.** Seven radial pieces with the breaks of `kernel.json`. A radius
   `r` in `(B_p, B_{p+1}]` belongs to piece `p`, and `r = 1` to piece 0. The
   support is the closed annulus `1 <= r <= 7/4`. Stage 1 checks every case on
   closed pieces, so it covers both one-sided values at each break
   (coordinator lines 237-247, `certificate_core.py` lines 104 and 107). My
   choice of the right-closed piece matches every count indicator
   `1[|x| <= R]` (closed balls) used by stage 1 (coordinator line 805,
   `signed_spherical_exact.py` line 85, `full_pair_count_gauge.py` line 42).
   On a closed cell the polynomial of piece `p` is the value of `K` on
   `(B_p, B_{p+1}]` and its continuous extension to `B_p`.
2. **Radial basis.** Piece `p` uses the degree-6 Bernstein basis on
   `radial_basis_intervals[p]` (`[1, 3/2]` for pieces 0 to 5, `[3/2, 7/4]` for
   piece 6), coordinates `7p..7p+6` of every 49-vector
   (`certificate_core.py` lines 41-53; `interval_certificate.py`
   `radial_bounds`).
3. **Gram origin data.** In `kernel.json` block 0 every ray has origin
   coordinate 0. The Gram origin entry is `kappa` and the origin row is `v`
   from `origin.json`, so `K_gram(x, 0) = v . phi(|x|)` and
   `K_gram(0, 0) = kappa`. Block 0 enters the nonorigin kernel without its
   origin row and column (coordinator lines 775-777 and 833;
   `certificate_core.py` line 69).
4. **Angular normalization.** `Z_l = U_l/(l+1)` (paper Section 3; checked
   independently: `certificate_core.py` `gegen` line 55 uses the closed
   formula, `kernel.py` the three-term recurrence).

## Triangle coordinates (the 112 cases)

5. **Full cases.** Coordinates are the side lengths `(a, b, c) = (|x|, |y|, |x-y|)`,
   regions an unordered triple `i <= j <= k`. `Gamma = k(a,b,t_c) + k(a,c,t_b) + k(b,c,t_a)`
   with `t_c = (a^2+b^2-c^2)/(2ab)` and so on. The certified polynomial is
   `(2abc)^18 Gamma_poly` (`certificate_core.py` lines 74-89,
   `signed_spherical_exact.py` line 170).
6. **Exterior cases.** Coordinates `(a, b, z)`, `z in [0, 1]`, with
   `c^2 = (a+b)^2 - z((a+b)^2 - (7/4)^2)`, so
   `t = -1 + z((a+b)^2 - 49/16)/(2ab)`; `z = 0` is `c = a+b` and `z = 1` is
   `c = 7/4` (`certificate_core.py` line 92, `signed_spherical_exact.py`
   line 135). The polynomial is `(2ab)^18 k(a, b, t)`. For `c > 7/4` the
   other two anchors fall outside the closed support, so `Gamma = k(a,b,t)`.
   At `z = 1` (`c = 7/4` exactly, inside the support) the true `Gamma` has
   two more terms; that point is a full case. The exterior polynomial is the
   limit `c -> 7/4+`.

## Which component enters which term

7. **Polynomial part of the cyclic sum.** Stage 1 builds
   `C/dd = (2abc)^18 (Gram + signed spherical) - eta (2abc)^18` with
   `eta = -(nested count + conditional count polynomial + conditional ray mixed term)`
   (coordinator lines 167-181; `certificate_core.py` line 111;
   `signed_spherical_exact.py` `build_signed` line 216). In the exterior case
   the third count flags are zero (coordinator line 170). `kernel.py`
   `gamma_poly` builds the same five components from their kernels:
   - signed spherical: `k = -beta g(t) chi(r) chi(s)`, `chi = 1[|x| <= 102597/100000]`,
     `h = beta B chi/2`, diagonal `-beta g(1) chi` (`origin.json`
     `signed_spherical_contract`);
   - count polynomial `q(N) = a + 2 d.N + N^T B N` with `N_i = sum 1[|x| <= R_i]`:
     `K(0,0) = a`, `h = d.n(x)`, `k = n(x)^T B n(y)` (coordinate order and
     symmetric upper triangle as in `base/new_work/spherical_sdp/count_quadratic.py`
     lines 1-9 and 50-62);
   - conditional ray `(23 - N_r)(N_R - N_r)`: `h = 23 beta 1[r < |x| <= R]/2`,
     `k = -beta (a(x) b(y) + b(x) a(y))/2` with `a = 1[|x| <= r]`,
     `b = 1[r < |x| <= R]` (COMPONENTS.md row 6; coordinator lines 247, 258, 790).
8. **Nonpolynomial corrections of the cyclic sum.** Caps, tents, allocations,
   the pair count gauge and the graph energies enter as the nonnegative
   amount `-sum k_c` subtracted from `Gamma` (lift helper
   `origin_route/interval_combined_binary_six_graphs_correlated_small_cap_bounds.py`
   line 140). Their kernels in `kernel.py`:
   - cap: `h = beta M w~(r)/2`, `k = -beta (w~(s) P(x|y) + w~(r) P(y|x))/2`,
     `P(x|y) = clip((1 - (c(r) - sigma t)/sqrt((1-sigma^2)(1-t^2)))/2)`,
     `c(r) = (r^2 + d^2 - a^2)/(2 d r)`, point mass at `sigma = +-1` and at
     `t = +-1` (`proofs/current/base/new_work/pointed/THEOREM.md` sections 2
     and 4). The pivot caps use `w~ = w 1[r > C]` (`origin.json`
     `cutoff_convention`: "actual pivot weight zero for r<=cutoff").
     Stage 1 uses the untruncated `w` in the cyclic sum
     (`base/new_work/pointed/interval_certificate_full_closed.py` lines 65-67)
     and in the origin row (coordinator line 781), and subtracts
     `J = beta w M` for `r <= C` in the pair term
     (`base/new_work/pointed/pair_interval_certificate_closed.py` line 46).
     The two kernels agree on feasible triangles by
     `TRUNCATED_PIVOT_THEOREM.md` section 2; the sampled comparison below
     tests this at every sample point;
   - tent: `h = beta C w(r)/2`, `k = -beta (w(s) phi(x + d yhat) + w(r) phi(y + d xhat))/2`,
     `phi(z) = (1 - 2|z|^2)_+`, `|x + d yhat|^2 = r^2 + 2 d r t + d^2`
     (`base/new_work/pointed/interval_certificate_tent_closed.py` lines 1-14
     and 73);
   - allocation: `h = beta u(r)/2`, `k = -beta (u(r) + u(s)) Cap(|x+y|^2)/2`,
     `|x+y|^2 = r^2 + s^2 + 2rst`, `u = w . phi` (`origin.json` `w`), two-level
     capacity `1` for `|x+y| < a`, `1/2` for `a <= |x+y| <= b`, else 0
     (`base/new_work/pointed/two_capacity_antipodal_allocation_exact.py`
     lines 44-48 and 116);
   - graph `j`: `K(0,0) = delta_j C_j`, `k = -delta_j 1[(x,y) selected]`, with
     closed radial and cosine windows from paper Table 1. Graph A uses an
     open left end at `103/100` for its far window
     (`six_exterior_graph_route/graph_repair_bounds.py` line 56); its weight is 0;
   - gauge: `h = w(r)/2`, `k = -(w(r) psi(s) + w(s) psi(r))/2`,
     `psi = 1[|x| <= 50379/50000]/24`, `w = -q0/(1 - psi)`
     (`base/new_work/pointed/full_pair_count_gauge.py` lines 1-18 and 142).
9. **Pair term.** Stage 1's pair polynomial `D` (coordinator line 853) is
   `2 h + k(x,x)` of every component except the gauge, with the caps entering
   through their untruncated linear term `beta M w/2`; the cap self term `J`
   is subtracted separately (`pair_interval_certificate_closed.py`
   `diagonal_bounds`). Tents, allocations, conditional rays and graphs have
   zero diagonal. `q0 = f + D - J`. The gauge then makes the full pair term
   `q0 + w(1 - psi) = 0` on `[1, 7/4]`. The comparison is of `q0`.
10. **Origin constant.** `A = kappa + nested a + sum beta_j a_j + sum delta_j C_j`
    (COMPONENTS.md "The objective uses these constants"; coordinator lines
    567 and 810). Caps, tents, allocations, spherical rays, conditional rays
    and the gauge have origin constant 0. `U_raw = (F(0) + A)/8`.
11. **Data fields not used.** The broad19 records in `origin.json` are not
    added: they are already Gram rays (COMPONENTS.md row 2). The
    `nested_count_quadratic.constant_rounding_repair` field is not added to
    the coefficients (the coordinator uses `coefficients` as stored, line 581).
    `functional_origin_v`, `pair_origin_v`, `objective` and all provenance
    fields are not read by `kernel.py`.
