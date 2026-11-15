# Component identities and boundaries

This note identifies the **same mathematical kernel** in `kernel.py` and the
stage 1 verifier. It complements exact coefficient comparisons; it does not
replace stage 1's interval sign proof. `contracts.py` checks the component
parameters taken from the live stage 1 worker state, not a saved report.
The nonpolynomial code still has a written-proof trust boundary. Its sampled
enclosure comparisons are diagnostics, never a proof on a continuum.

All paths below start at `../stage1_v7/record/current/`. Write `r=|x|`,
`s=|y|`, `t=<x,y>/(rs)`, and `w(r)` for the physical Bernstein polynomial
on the selected radial piece. Stage 1 uses the closed extension of every
radial piece. The actual kernel selects `(B_p,B_(p+1)]`, with `r=1` in piece
zero; thus every actual value is among the certified closed extensions.

## Polynomial part and origin

The 49 Bernstein basis polynomials are compared by exact coefficients with
`certificate_core.radial`. The independent code forms Gram matrices by
FLINT integer matrix multiplication, and generates `U_l/(l+1)` by its
three-term recurrence. It then expands Gram, signed-spherical, nested-count,
conditional-count and conditional-ray terms from their individual kernels.
For each of the 84 full and 28 exterior cases every coefficient of the
homogenized sum is compared with stage 1's actual `build_signed` output.
The factors are `(2abc)^18` and `(2ab)^18`, respectively, and are positive
on the entire domain. Coefficient identity therefore proves identity of
these polynomial parts on every real point of the branch, including its
boundary. The independent program does not import stage 1's assembler.

The origin is

`A = kappa + nested_constant + sum(weight * conditional_constant)
     + sum(graph_weight * graph_capacity)`.

All other components have zero origin constant. The checks compare `A`,
`(F(0)+A)/8` and the exact census objective as rationals, and independently
verify `M_0 u=v`, `u.v=kappa`. They do not treat the data's reported objective
as a derivation.

## Caps, including pivot truncation

Set `c(r)=(r²+d²-a²)/(2dr)`, `A=c(r)-sigma*t` and
`S=(1-sigma²)(1-t²)`. For the conditional direction on a two-sphere the
remaining coordinate is uniform on `[-1,1]`. Hence the closed-ball event
has probability

* `1[A<=0]` if `S=0`;
* `clip((1-A/sqrt(S))/2,0,1)` if `S>0`.

This is exactly the function enclosed by
`base/new_work/pointed/interval_certificate_general_closed.py:probability_bounds`.
The `S=0,A=0` value is **one**, including `t=±1` or `sigma=±1`; the ratio
formula is not used there. `kernel.py:cap_probability` uses this case split.
Squaring is used only with the sign of `A` known, so its tests for zero and
one do not introduce extraneous solutions.

With `wtilde(r)=w(r)1[r>C]` for a truncated cap (and `wtilde=w` otherwise),

`h(r)=beta M wtilde(r)/2`,
`-k(x,y)=beta (wtilde(s)P(r,t)+wtilde(r)P(s,t))/2`.

`contracts.py` compares `d,a²,sigma,M,beta,w,C` for every cap, and checks
the strict onset inequality used by stage 1. In that case `sigma=1`,
`d²<=a²<1/2`, `C>=1` and `d C < 1+d²-a²`. If `1<=s<=C` and `x,y` are
distinct admissible points, then

`|x-d yhat|² >= (1-d/s)|x|² + d²-ds+d/s >= 1+d²-ds > a²`.

Here `d<1<=s` justifies the second inequality. Therefore the probability
is zero whenever truncation changes a pivot weight. The untruncated cap
correction used by stage 1 equals the actual truncated one on every
admissible triangle, including `s=C`.

On the diagonal, `P(r,1)=1[r²-2d sigma r+d²-a²<=0]` exactly. Let `J` be
the reduction stage 1 subtracts from its untruncated linear pair term.
For each cap,

`J_cap = beta w(r) M` if `r<=C`, and
`J_cap = beta w(r) 1[r²-2d sigma r+d²-a²<=0]` otherwise.

For an untruncated cap only the second expression applies. Subtracting
this `J_cap` from `beta M w(r)` gives `2h(r)+k(x,x)` in both cases.
This proves the pair identity at every radius, including irrational roots
of the displayed quadratic and the pivot cutoff. The equality uses the
same parameters checked by `contracts.py`; the extra `J` samples in
`pair.json` are not the proof and do not cover irrational endpoints by
themselves.

## Tents

The tent is `phi(x+d yhat)=(1-2(r²+d²+2drt))_+` and

`-k(x,y)=beta (w(s)phi(x+d yhat)+w(r)phi(y+d xhat))/2`.

The quantities `d,beta,w,capacity` agree exactly with the validated stage 1
records. Substituting `2rst=r²+s²-c²` gives the expression in
`interval_certificate_tent_closed.py:full_tent_bounds`; substituting
`t=-1+z((r+s)²-(7/4)²)/(2rs)` gives its `outside_tent_bounds` expression.
The hinge is continuous at zero, so no endpoint choice remains. Stage 1's
full-branch shortcut is valid: positive tent weight forces
`c²>2+sqrt(2)>(7/4)²`. Equivalently one can keep the zero hinge explicitly.

## Allocations

Both kernels use `-k(x,y)=beta(u(r)+u(s))Phi(|x+y|²)/2`.
The validated stage 1 profile is `u(r)=(2r-2)^4` on its inner basis and
zero on the outer basis. Equality of its Bernstein coefficients is checked
exactly, including both one-sided values at `r=3/2`.

For the classical profile, the function is

`Phi(q)=min(1,2(1-sqrt(q))²)` for `0<=q<1/4`,
`Phi(q)=1-2q` for `1/4<=q<1/2`, and zero for `q>=1/2`.

The pieces agree at both joining points. For the two-level profile the
value is one for `q<a²`, one half for `a²<=q<=b²`, zero for `q>b²`.
The exact `a,b` and the strict/closed conventions are checked against the
stage 1 validated record. `|x+y|²=r²+s²+2rst`; the exterior substitution
is the squared-distance expression used by the stage 1 helper. On full
branches `|x+y|²=2r²+2s²-c² >= 4-(7/4)²>1/2`, so both profiles vanish.
The allocation boundary checks include `a²,b²,1/4,1/2` themselves.

## Graphs

Every graph component is `-k=delta 1[radial windows and angular window]`.
The six weights are compared exactly with stage 1. Graphs A and B have
weight zero. The four active graphs have closed radial and angular
windows; the radial and angular interval endpoints and intervening cells
are exercised by exact rational singleton checks. Every selected pair has
distance strictly larger than `7/4`: this follows by evaluating
`r_min²+s_min²-2r_min s_min t_max` for each active graph (both radii increase
this expression because `t_max<0`). This exact inequality is checked.
Thus these graph terms vanish on full branches, as stage 1 assumes.

The graph checks test the actual stage 1 exterior coordinate transformation
and all relevant indicator endpoints. Their interpretation as an exhaustive
indicator check uses the written definitions in
`six_exterior_graph_route/{graph_repair,outer_exterior_graph,narrow_exterior_graph,moderate_exterior_graph}_bounds.py`.

## Pair gauge

Write `q0(r)=exp(-alpha*r²)F(alpha*r²)+D(r)-J(r)` for the ungauged pair
term. Exact coefficient checks identify `F,D`; the cap identity above
identifies `J` on all radii. The code compares the scalar exponent, gauge
radius, count bound, all seven gauge flags and the absence of additional
cuts with the actual `FullPairCountBudget` used by stage 1.
Therefore both implementations use exactly

`psi(r)=1[r<=R]/M`, `w_g(r)=-q0(r)/(1-psi(r))`,
`-k_g(x,y)=(w_g(r)psi(s)+w_g(s)psi(r))/2`.

Since `M>=2`, the denominator is positive. The stage 1 fresh pair proofs
establish `q0<=0` before its nonnegative clamp is enabled. That clamp only
tightens an enclosure and changes no mathematical value. The resulting
full pair term is identically `q0+w_g(1-psi)=0`, including the closed gauge
radius and every jump of `J`.

## Extent of the conclusion

The exact polynomial comparisons plus these component identities identify
the independent definition with the kernel certified by stage 1. The
coefficient equalities, parameters and specified endpoint tests are machine
checked. The interpretation of the interval routines as enclosures of the
displayed cap/tent/allocation/gauge functions remains a written argument
and code-review obligation, not a Lean theorem or a second independent
continuum sign proof. Stage 5 has its separate symbolic link check. No
sample agreement is promoted to a universal identity.
