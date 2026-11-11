# Dimension-sensitive lens caps

This extends the elementary lens proof without changing any frozen source or
running a density optimization. All code certificates and tail budgets below
have exact local replay checks. They are not themselves density bounds.

## Strictly obtuse vectors

If vectors v_1,...,v_M in R^m satisfy v_i dot v_j<0 for every i!=j, then
M<=m+1. Indeed, if M>m+1 they have a nonzero affine dependence

    sum_i c_i v_i=0, sum_i c_i=0.

There are both positive and negative coefficients. Put
w=sum_(c_i>0)c_i v_i=sum_(c_j<0)(-c_j)v_j. Then

    |w|^2=sum_(c_i>0,c_j<0)c_i(-c_j)(v_i dot v_j)<0,

a contradiction. This proves the dimension bound without assuming vector
independence, unit lengths, a positive dependence, or a particular rank.

For common neighbors centered at the base midpoint, |v_i|^2<=rho^2 and
separation at least 1 imply v_i dot v_j<=rho^2-1/2. When rho^2<1/2, the cap is
n+1. At c=41/25 and r=27/25, rho^2=247/500, so the dot upper is -3/500.
It decreases for larger c. Taking the minimum with each centroid cap yields
`step_allowances_midpoint_rankin_005.json`.

## A stronger transverse cap

Write a common neighbor as (ell,z) relative to the base midpoint and direction.
Both endpoint squared distances lie in [1,r^2], so

    |ell| <= (r^2-1)/(2c).

Let A=(r^2-1)^2 and D=A/c^2. For two common neighbors, the longitudinal
difference squared is at most D, and hence their transverse separation squared
is at least 1-D. Also |z|^2<=U=r^2-c^2/4. Therefore

    z_i dot z_j <= U-(1-D)/2.

At c=41/25 this upper bound is exactly -3583/4202500<0. Both U and D decrease
with c, so it stays negative for all larger c. The z vectors lie in R^(n-1),
so the strictly obtuse argument gives M<=n. No normalization is needed. A zero
transverse vector is harmless: it cannot coexist with another vector under
the strict negative off-diagonal dot bound.

This gives the strictly stronger artifact
`step_allowances_transverse_rankin_005.json`. Its early caps are 5,6,7 rather
than 6,7,8 in dimensions 5,6,7. Both artifacts preserve the old centroid
artifact separately and all 312 budgets in each pass Arb replay.

## Spherical codes for c below 1.64

### Midpoint projection

Every common neighbor v satisfies |v|<=rho<1 for the base range here. If v=0
occurs, no second neighbor is possible: its distance from v would be at most
rho<1. Otherwise normalize every v to the unit sphere. For positive lengths
s,t<=rho and separation at least 1, their normalized dot is at most

    (s^2+t^2-1)/(2st) <= 1-1/(2rho^2).

The derivative with respect to s has sign s^2-t^2+1, which is at least
1-rho^2>0; the other derivative is analogous. This proves the inequality
without a positive lower norm assumption. Projection is injective: two vectors
on the same positive ray have distance less than 1. The angular upper bound
decreases with c, so a code bound proved at a base threshold remains valid for
larger c. The exceptional zero case already has M<=1.

### Transverse projection

For the transverse vectors above, squared norms lie in

    L(c)=1-c^2/4-A/(4c^2) <= |z|^2 <= U(c)=r^2-c^2/4.

Throughout c in [3/2,41/25], L(c)>0. It decreases on this range, so checking
the positive value at 41/25 suffices. The separation bound gives normalized
dot at most

    (u+v+D-1)/(2sqrt(uv)),  u,v in [L,U].

This increases in each argument whenever L-U+1-D>0. Here that margin is
2-r^2-5A/(4c^2), which increases with c and is positive already at 3/2.
Thus the angular upper is

    s(c)=1-(1-A/c^2)/(2(r^2-c^2/4)).

This decreases with c. To verify the last assertion algebraically, write
x=c^2 and q=r^2. The derivative of the subtracted quotient has positive
numerator

    x^2-2Ax+4qA = (x-A)^2 + A(4q-A)>0.

All denominators are positive in this range. Consequently a transverse code
certificate at a base threshold is valid until 41/25, where the stronger
obtuse cap takes over. Norm positivity also excludes the zero transverse
vector; the angular bound s<1 makes normalization injective.

The code records include exact norm, angular, monotonicity, and uniform-domain
quantities. Replay reconstructs these, including the positive lower transverse
norm at 41/25 and the positive derivative remainder A(4q-A).

### Exact Delsarte certificates

For a code in R^m with off-diagonal inner products <=s, the stored polynomial
Q(t)=sum_l a_l P_l^(m)(t) has a_0=1, a_l>=0, and Q<0 on [-1,s]. The addition
theorem yields M^2<=sum_(i,j)Q(t_ij)<=M Q(1). Therefore
M<=floor(Q(1)). The replay checks all coefficient signs, reconstructs the full
polynomial exactly, checks negative endpoint signs, uses exact rational Sturm
root counting to prove negativity on the entire interval, and reconstructs
the rational objective and its integer floor. Stored success booleans are
never used as evidence.

The exact results at c=1.50 / 1.55 / 1.60 are:

| Packing dimension | Midpoint caps | Transverse caps |
|---|---|---|
| 5 | 14 / 12 / 11 | 10 / 9 / 8 |
| 6 | 19 / 16 / 13 | 15 / 13 / 11 |
| 7 | 24 / 20 / 16 | 20 / 16 / 13 |

The transverse bounds are stronger in all nine cases. The 24 preserved code
inputs also include both projections at 1.64, but the transverse obtuse cap
5/6/7 is stronger there. The inputs and replay are `spherical_lens_caps.json`,
`spherical_lens_caps.py`, and `spherical_lens_caps_replay.json`.

## Combined portable allowance artifact

`step_allowances_spherical_005.json` combines the transverse spherical caps
from 1.50 to 1.64, transverse obtuse caps thereafter, and sharper centroid caps
when available. It has 132 width-.005 bins per dimension through 2.16. Each
bin allowance is approximately 85% of the available physical-tail budget and
is verified independently by Arb. All 396 bins passed with no fallback
reductions. Replay checked 162/182/168 Taylor cells in dimensions 5/6/7, with
maximum depths 6/7/5.

The artifact includes every spherical polynomial needed for its three early
thresholds in each dimension. Its verifier reads no external certificate or
original source: it reconstructs the geometry, spherical certificates, exact
physical coefficient payload, partition, boundary minima, and analytic pair
budgets from the artifact. Source-wide SHA values are provenance only.

At exactly c=R=3/2 the allowance is ZERO, preserving the original compact
support boundary. On the open first bin c>R only the far kernel survives.
At every other bin endpoint use the smaller neighboring allowance; the
existing exact interval-query helper implements this convention.

`construct_step_allowances.py` now exports `validate_cap_geometry(record)` for
the far verifier to call instead of a centroid-only test. It verifies embedded
polynomials and exact formulas, and returns the validated cap table. Full
`--verify` additionally proves every budget and checks the physical hash.

Reproduction:

    python construct_step_allowances.py --geometry spherical_transverse --output step_allowances_spherical_005.json
    python construct_step_allowances.py --verify step_allowances_spherical_005.json --candidate PATH_TO_FIXED_F_CANDIDATE

For portable replay copy `construct_step_allowances.py`, `spherical_lens_caps.py`,
the budget JSON, and optionally `step_allowance_query.py` together. The replay
requires Python with python-flint only. NumPy/SciPy are imported lazily for
separate spherical-code construction and are not needed for replay. Replay
does not invoke a numerical LP solver. The frozen original budget artifact
still replays.
