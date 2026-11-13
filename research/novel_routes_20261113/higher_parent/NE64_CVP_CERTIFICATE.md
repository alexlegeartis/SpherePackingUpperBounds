# Explicit Ne64 dual coset and certified 65-dimensional lattice

## Status

Class 4 of the deterministic search has an exact vector of squared norm `26/3` in `Ne64*`, and both complete discovery traversals found no vector of squared norm below `26/3` in that coset. The coefficient-logged run visited exactly **4,567,808,611 nodes**, finished without its 20-billion-node safeguard, and returned `empty_float_candidate` in 179.251 seconds. The name of that raw status is intentionally retained.

`audit_ne64_rounding.py` proves, with rational arithmetic, a uniform error margin for the logged floating traversal. The argument connecting that audit to exhaustive enumeration is given below. The independent outward-interval enumeration has now completed in all four shards, covering 4,880,529,145 nodes with zero leaves. Its summary is `../ne64_hole_certificate/summary.json`, and it is the primary certificate. The generated 65-dimensional lattice is certified; `verify_ne65.py` checks its algebra and the certificate inputs. This document does not assert world-record priority.

## Exact arithmetic and the certified lattice

The published Ne64 lattice `L` is even, has minimum 12 and determinant `3^32`, and is 3-modular. The archived original Gram is `G`. The input file `ne64_primal_best.npz` supplies exact integer matrices `G,H,T,Ti,R` satisfying

\[
GH=3I,\qquad TTi=TiT=I,\qquad R=TGT^T.
\]

The class is represented in reduced coordinates by the residue in `coset_full_4.txt`; its original and reduced residues are also in `coset_classes.json`. If its reduced representative is `c`, enumeration asks whether

\[
z\equiv c\pmod3,\qquad zRz^T\le60
\]

has a solution. The original numerator `w=zT` belongs to the specified dual class because `wG` is divisible by 3. Every numerator in the class has norm congruent to6 modulo18. Thus the next possible value after 60 is78; absence through 60 means the physical coset minimum is at least `78/9=26/3`.

`coset_class4_norm78.jsonl` contains a witness at exactly 78, found in 1,305,912 nodes. The witness was independently checked in both bases and reduced modulo 3. Its original numerator `w` gives the rational dual vector `h=w/3`. Set

\[
t=\sqrt{10/3},\qquad L_{65}=\langle (L,0),(h,t)\rangle.
\]

Its saved Gram is

\[
G_{65}=\begin{pmatrix}
G & Gw^T/3\\
wG/3 & 12
\end{pmatrix}.
\]

Every entry is integral and the diagonal is even. The Schur complement is `10/3`; exact determinant evaluation gives

\[
\det G_{65}=10\cdot3^{31}=6176733962839470.
\]

The certified coset bound gives minimum 12. In layers not divisible by 3, the horizontal squared length is at least `26/3` and the vertical contribution is at least `10/3`. Nonzero layers divisible by 3 have vertical squared length at least 30. Layer zero inherits the published minimum 12.

Consequently the center density is

\[
\delta_{65}=\frac{3^{17}}{\sqrt{10}}
=40,837,705.2485403\ldots.
\]

This exceeds the earlier existential depth8 extension by the factor `sqrt(6/5)`, or 9.5445115 percent. The compact construction is in `ne65_candidate.npz`, its integer Gram in `ne65_candidate_gram.txt`, and exact invariants in `ne65_candidate.json`.

## Logged floating traversal and error model

The executed source is `coset_enum_logged.cpp`, compiled with

```sh
c++ -O3 -std=c++17 -fno-fast-math coset_enum_logged.cpp -o coset_enum_logged
```

The implementation uses ordinary IEEE754 binary64 arithmetic, round to nearest. Long-double initialization does not enter the proof: every stored binary64 coefficient used in enumeration is logged in hexadecimal, and compared with the exact rational LDL decomposition. Unsafe floating reassociation is not enabled. Fused operations, if generated, have no more error than the conservative bounds below.

Write

\[
R=LDL^T,\quad U=L^T,\quad
zRz^T=\sum_i d_i\left(z_i+\sum_{j>i}U_{ij}z_j\right)^2.
\]

Let `u=2^-53` and `gamma(k)=ku/(1-ku)`. The audit parses the hexadecimal values as exact dyadic rationals. All displayed decimal bounds below are explanatory; assertions use rational numbers.

For any true target vector of norm at most 60, Cauchy–Schwarz gives

\[
|z_j|\le M_j:=\left\lceil\sqrt{60(R^{-1})_{jj}}\right\rceil,
\qquad
\left|z_i+\sum_{j>i}U_{ij}z_j\right|
\le Z_i:=\left\lceil\sqrt{60/d_i}\right\rceil.
\]

The inverse is obtained exactly as `R^-1=Ti.T H Ti/3`. The largest `M_j` is 416.

A center uses four accumulators. Fewer than 70 rounded operations occur along the path of any summand; the standard absolute dot-product bound is therefore conservatively `gamma(70)` times the sum of absolute products. If `uf` denotes the stored coefficients, the exact bound is

\[
\eta_i=\sum_{j>i}|uf_{ij}-U_{ij}|M_j
+\gamma(70)\sum_{j>i}|uf_{ij}|M_j.
\]

It includes both coefficient error and all center arithmetic. Its maximum is less than `1.533045e-11`.

The program's squared radius is the binary64 constant

\[
R_0=60+\frac{8796093}{8796093022208},
\]

so the exact excess over 60 is about `9.999999974752427e-7`.

## Closest-base rounding, ordinary branches, and accumulated sums

The nearest allowed coordinate is selected by

```cpp
mid=nearbyint((center-c[i])/3.0);
base=c[i]+3*mid;
delta=base-center;
```

The quotient has two rounded operations. If the computed center has absolute bound `Chat_i`, define `eps_i=gamma(2)*(Chat_i+2)`. The real distance between `base` and the computed center is at most `1.5+eps_i`. Rounding the quotient can select the adjacent integer at a near tie, but the excess squared distance over the exact nearest allowed integer is at most

\[
6\,eps_i+4\,eps_i^2.
\]

This near-tie penalty is explicitly included in the audit; exact-nearest selection is not assumed.

For the intended target coordinate, a step is formed by adding an exact multiple of3 to the rounded `delta`. The audit bounds the extra subtraction and addition error by

\[
\kappa_i=u(1.5+eps_i)+u\bigl(Z_i+\eta_i+u(1.5+eps_i)\bigr).
\]

It then bounds each computed squared term using the exact stored diagonal error, `eta_i+kappa_i`, and `gamma(2)` for the two multiplications. It sums these errors and includes `gamma(64)` for accumulated nonnegative partial sums. The separate initial-base bound includes the near-tie penalty above.

At a node on any true target path, the remaining-radius subtraction has absolute error at most `61u`. The resulting rational bounds are:

| Test | Uniform possible overestimate |
|---|---:|
| Target-coordinate branch and accumulated partial sum | `1.012023691112e-8` |
| Initial nearest-base test and accumulated partial sum | `1.015165711089e-8` |

Both are strictly smaller than the exact radius excess. Hence a true norm-at-most60 path cannot be rejected by the remaining-radius check, the initial-base check, or its own branch test. The residual safety margin exceeds `9.8984834e-7`.

## Why the early stopping rule does not skip a target branch

The loop tries `delta+3*step*sign` and `delta-3*step*sign` for integer `step>=1`, and stops when neither is admitted. The audit proves `abs(delta)<1.6`. On each side, from step 1 onward, the exact absolute values strictly increase. IEEE rounding is monotone, as are multiplication by a positive diagonal and the resulting squared magnitude. Therefore the tested squared terms on each side are nondecreasing.

If a target coordinate occurs at step `s`, its own branch passes by the preceding error bound. Every earlier step on that same side also passes. The `if(!one) break` condition therefore cannot terminate the loop before the target step. The initial base was separately shown to pass, including near ties. This proves branch completeness by induction from coordinate 63 down to 0.

The full traversal finishes with no node cap. At every reached leaf, the numerator norm is evaluated with exact signed128-bit integer arithmetic. No leaf with norm at most 60 was returned.

## Arithmetic-domain checks

The audit separately bounds all visited coordinates, without assuming they are true norm-at-most60 vectors. A downward triangular recurrence using the stored coefficients, the positive stored diagonals and radius61 gives a bound of449,695,303. The extra unit in that recurrence dominates all center and step rounding errors; the audit explicitly checks the global roundoff terms are below one tenth.

These bounds keep integer-to-double conversions exact, the nearby-integer result and coordinate arithmetic inside signed64-bit range, and the exact leaf quadratic form inside signed128-bit range. Exact dyadic denominator bounds keep all nonzero arithmetic intermediates above the subnormal range. The maximum relevant dyadic denominator exponent is recorded in `ne64_rounding_audit.json`.

Thus the pruning argument does not assume exact floating arithmetic, exact quotient ties, or an unchecked absence of integer overflow.

## Reproduction and audit artifacts

```sh
# Execute the complete logged traversal.
./coset_enum_logged ne64_primal_best_Q.txt coset_full_4.txt 20000000000 ne64_certified_run_ldl.hex

# Check exact coefficient errors and all uniform bounds.
python audit_ne64_rounding.py
```

The audit stores SHA256 hashes of the matrix, input class, logged coefficients and executed source. `coset_full_4_logged.jsonl` is the recorded complete-run result. The separate capped exploratory screens are excluded from the checkpoint. Their `node_cap` entries remain unresolved and are not part of the certificate.
