# A Gram-assisted far-pair energy refinement in dimension four

The result is **upper packing density < 0.629520189662** for every packing
of congruent balls in R4. It improves the strict bound 0.629521681653 in
[Felpix-Studios/sphere-packing](https://github.com/Felpix-Studios/sphere-packing/tree/bb97550ba47ad92bacfee567843888be16895a23/preprints/D4).
This is an extension of that certificate. Its scalar function, pair kernel,
count census, graph terms and density-transfer theorem are retained.

All upstream files in `vendor/D4/` are unchanged copies at commit
`bb97550ba47ad92bacfee567843888be16895a23`. Their authors retain credit for
the underlying 0.629521681653 proof; see [NOTICE.md](NOTICE.md).

## 1. An unused positive energy

Use the upstream normalization: minimum centre distance 1, ball volume
pi^2/32, Fourier integral pi^2/4 and support radius R=7/4. Let
Z_d(t)=U_d(t)/(d+1) be the normalized S3 Gegenbauer polynomial. For each
positive angular degree, the original energy contains

    E_d = sum_(x,y in X) Z_d(xhat.yhat) phi(|x|)^T M_d phi(|y|),

where M_d is the sum of positive integer Gram rays divided by their common
positive denominator. These are row 1 of the upstream
[component inventory](vendor/D4/verify/stage1_v7/COMPONENTS.md).

The census proof uses the degree-zero block, count energies and two
separately appended spherical energies. It drops every E_d with d>=1 as
nonnegative; see its [projection and minimization proof](vendor/D4/verify/stage2_origin/PROOF.md)
and [finite-count refinement](vendor/D4/verify/stage3_census/finite_count_gain.md).
Thus the same census floor remains valid after retaining explicit lower
bounds on selected E_d. The broad19 rays are already inside M_d and are
counted once. No auxiliary spherical proof is reused as an extra energy.

Put R1=501473/500000. In the 49-coordinate radial Bernstein basis, let v
have seven initial entries 1 and every other entry 0. The closed-inner
radial convention gives v.phi(r)=1 exactly when 1<=r<=R1, and 0 otherwise.
For d=5,6,7,8,9 the exact records in `gram_credits.json` give

    M_d w_d = v,   b_d = v^T w_d > 0,   beta_d = 1/b_d.

Consequently M_d-beta_d*v*v^T is positive semidefinite. For example, with
T=I-beta_d*w_d*v^T, it equals T^T*M_d*T. Schoenberg positivity then gives

    E_d >= beta_d S_d,
    S_d = sum_(i,k=1..N) Z_d(u_i.u_k) >= 0,

for the N near neighbours within R1, even in the presence of additional
neighbours elsewhere in the support. `gram.py` independently rebuilds M_d
from the original rays and checks all 245 range equations exactly. The
approximate beta values are 6.9553e-9, 3.1421e-9, 2.5660e-10, 1.0141e-9
and 4.8598e-11. The rational records, rather than these displays, are used.

## 2. A polynomial that shares those Gram credits

Let a=133/250 and let p(t)=sum_(d=0..16) c_d Z_d(t), including its scale,
be the rational polynomial in `polynomials.json`, entry `1`.
The new checks establish

    p(t) <= 0                 for -a <= t <= s1,
    p(t) <= -K(y,z)           for |y|,|z| in [1,R1], t in [-1,-a],
    c_d + beta_d >= 0         for every d>=1,

where s1=1-1/(2*R1^2), and beta_d=0 outside the five degrees above.
The c_d at degrees 5 through 9 are negative. The argument therefore uses
the retained Gram energies together with the reallocation energy; it does
not assert that p alone has positive Gegenbauer coefficients.

`polynomial.py` converts the normalized Gegenbauer recurrence to monomials
with Fractions and proves the central sign by exact Bernstein subdivision
on the whole closed interval. Ten leaves suffice. `boxes.py` replays all
308 leaves of a midpoint subdivision of the complete three-variable far
domain. It uses the unchanged outward-Arb lower-envelope engine and kernel
formula. Coverage, endpoint flags and every leaf are checked afresh; saved
producer minima are ignored. The upstream full replay separately checks
the exact link of this kernel formula to its cyclic sign certificate.

## 3. The joint lower bound

The upstream reallocation H_R is nonnegative. For near neighbours, unit
separation gives t_ik<=s1. Pairs with t_ik<-a have mutual distance greater
than 7/4, since 2+2a>49/16. At such a pair its reallocation term is exactly
-K(y_i,y_k). Dropping all other nonnegative terms gives

    H_R >= sum_(t_ik<-a) p(t_ik) >= sum_(i!=k) p(t_ik),

because p<=0 on every remaining allowed angle. With S_d as above,

    H_R + sum_d beta_d S_d
      >= N^2 c_0 - N p(1) + sum_(d>=1) (c_d+beta_d) S_d
      >= N^2 c_0 - N p(1).

The left side is also nonnegative. Therefore the original census floor
for a group (N,1) can be increased by

    max(N^2 c_0 - N p(1), 0).

Nothing in K, the pair signs or the cyclic signs has changed. In particular,
the Gram credits are not subtracted from the kernel used for the far-box
check. They are lower bounds on energy that the previous census discarded.
This avoids consuming the same credit twice.

For groups with first containing radius j=2,3,4,5, retain the original
certified quartic and coefficient lambda5=1/1000000. Its contribution is
lambda5*max(N*(13561*N-179244)/250000,0). The inherited B5 certificate is
replayed in the full upstream proof. No larger-radius experimental
polynomial is used in this result.

## 4. Exact composition and density

Minimize the raised floors over exactly the original 124 groups, including
the empty neighbourhood. The binding group is (24,1), with increment

    0.000054047284957656

and resulting floor approximately 0.003837814256786624. The next group is
(26,5), at approximately 0.003838744001707000. Keep the original broad-graph
saving 3/8000000 in density units, once. The upstream finite-to-density
argument yields

    upper density <= U_raw - min_group_floor/8 - 3/8000000
                  < 0.629520189662
                  < 0.629521681653.

The unrounded objective is approximately 0.629520189661472. Its exact
rational value and input hashes are in `certificate.json`; `verify.py`
recomputes them from the selected census and the rational polynomial.
The reduction between the two strict decimal bounds is exactly
1491991/10^12. This establishes the stated comparison with the pinned
repository, without asserting an exhaustive worldwide-priority review.

## Trust and replay scope

The full verifier replays the unchanged upstream seven-stage chain,
including all 112 cyclic branches, the complete count domain, spherical
proofs, graph proof, both old box certificates and both kernel identity
checks. It then checks the new polynomial, exact Gram identities, all new
box leaves and the 124-group ledger. Written mathematical inputs remain
the upstream local-energy and density-transfer arguments, the addition
theorem, interval lower-envelope soundness, and the elementary Schur
identity proved above. There is no proof-assistant formalization claim.

The default command performs a fresh full replay. `--base-replay DIR`
explicitly reuses an existing full replay and records that fact.
`--extension-only` certifies only this extension against committed base
inputs and reports `PASS_EXTENSION_ONLY`; it is not a fresh global replay.
