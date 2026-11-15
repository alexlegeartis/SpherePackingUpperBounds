# Classical Blichfeldt improvement and movable antipodal centers

This is a local-energy construction and a finite-dual search diagnostic.
It is not a density bound. The radial capacity function is classical:
equation (5) in G. Fejes Toth, F. Fodor, and V. Vigh, *The packing density
of the n-dimensional cross-polytope*, gives it at separation two. The
formula below is rescaled to separation one.

Source: https://www.math.u-szeged.hu/~vigvik/preprints/blichfeldt.pdf

## 1. Capacity, with an elementary proof

Write `P(r)=(1-2r^2)_+`. Define

\[
P_*(r)=
\begin{cases}
1-P(1-r)=\min\{1,2(1-r)^2\},&0\le r\le1/2,\\
P(r),&r\ge1/2.
\end{cases}
\]

For every finite unit-separated set in any Euclidean dimension,

\[
\sum_i P_*(|z_i|)\le1.
\]

If every norm is at least `1/2`, this is the usual Blichfeldt inequality.
For completeness, on the `n` positive summands that inequality follows
from

\[
n\sum_i|z_i|^2-\left|\sum_i z_i\right|^2
=\sum_{i<j}|z_i-z_j|^2\ge n(n-1)/2.
\]

There is at most one point with norm below `1/2`. Let it be `z_0`,
with norm `a`. If there are no other positive ordinary summands, the
claim is immediate. Otherwise retain those other `m>=1` summands and put

\[
T=\sum_{j=1}^m|z_j|^2,\qquad t=\left|\sum_{j=1}^m z_j\right|.
\]

Their mutual separation and separation from `z_0` give, respectively,

\[
T\ge\frac{m-1}{2}+\frac{t^2}{m},\qquad
T\ge m(1-a^2)-2at.
\]

These imply

\[
T\ge\frac{m-1}{2}+(1-a)^2.\tag{1}
\]

Indeed the first bound suffices if `t>=sqrt(m)(1-a)`. Otherwise use the
second bound at this upper endpoint for `t`; its excess over the
right side of (1) is

\[
(m-1)(1/2-a)+a(1-a)(\sqrt m-1)^2\ge0.
\]

Thus the sum of the other positive ordinary terms is at most
`1-2(1-a)^2`. If that expression is negative, no such other positive
terms exist. Their contribution is therefore at most `P(1-a)`, which
is exactly the remaining capacity after the contribution `P_*(a)`.
All remaining points have zero contribution. This proves the claim,
including the endpoint cases by the displayed definition.

The improvement is pointwise: for `r<=1/2`, its excess over `P(r)` is
`min(2r^2,(1-2r)^2)>=0`. Its plateau is

\[
P_*(r)=1\quad\text{when }0\le r\le a_0:=1-1/\sqrt2.
\]

## 2. The reflected test point can move

Let `C` be any actual packing neighborhood in the annulus
`1<=|x|<=R`, with `R=7/4`. Let `u>=0` be radial, and choose a fixed
rational number `delta` satisfying `|delta|<=29/100`. Define

\[
c_\delta(y)=-y+\delta\frac{y}{|y|},\qquad
E_{u,\delta}(C)=\sum_{y\in C}u(|y|)
 \left[1-\sum_{x\in C}P_*(|x-c_\delta(y)|)\right].
\]

Translation preserves unit separation, so each bracket is nonnegative.
This energy has zero origin constant, origin coupling `u/2`, and
symmetric neighbor kernel

\[
K_{u,\delta}(x,y)=-\frac12\left[
u(|y|)P_*(|x-c_\delta(y)|)+
u(|x|)P_*(|y-c_\delta(x)|)\right].
\]

Both nonzero arguments must lie in the active annulus; otherwise this
kernel is zero. Its diagonal vanishes because the self distance is
`2|y|-delta>=171/100>1/sqrt(2)`. Thus its complete unlifted pair cost is
exactly `u(r)`. Adding `gamma E`, with `gamma>=0`, must retain the pair
constraint `q+gamma u<=0`, and any subsequent basic count gauge must
use that complete modified pair expression.

If `C` is centrally symmetric, the target `x=-y` has distance `|delta|`
from the test center and contributes one. Every other target has
distance at least `1-|delta|>=71/100>1/sqrt(2)` and contributes zero.
The exact squared comparison is `(71/100)^2=5041/10000>1/2`.
Consequently `E_{u,delta}(C)=0` on every centrally symmetric neighborhood,
including every radial truncation of every lattice packing. No
assumption that an arbitrary maximizing packing is centrally symmetric
is made or needed.

## 3. Only exterior triangles receive this negative kernel

For a full triangle let `r=|y|>=1`, `s=|x|>=1`, and `|x-y|<=R`. Set
`t=1-delta/r>0`. Then

\[
|x+t y|^2=(1+t)s^2+t(1+t)r^2-t|x-y|^2
\ge(1+t)^2-tR^2.
\]

At `R=7/4`, minimizing over all `t>=0` gives

\[
1-\frac{(R^2-2)^2}{4}=\frac{735}{1024}>\frac12.
\]

Hence every directed capacity term vanishes on every full triangle,
including all closed branch boundaries. The resulting full-triangle
row change is solely the basic count-budget cost already described in
`INDEPENDENT_ANTIPODAL_ALLOCATION.md`.

On an exterior branch use the existing parameter `0<=z<=1` with

\[
\cos\theta=-1+z\frac{(r+s)^2-R^2}{2rs}.
\]

For the pivot of radius `r`, the exact directed squared distance is

\[
|x-c_\delta(y)|^2=(s-r+\delta)^2+
 z\left(1-\frac\delta r\right)((r+s)^2-R^2).
\]

The swapped pivot has its own formula with `r` and `s` exchanged.
They must not be replaced by one common capacity value when `delta!=0`.
All denominators are separated from zero by `r,s>=1`.

## 4. What the finite-dual test found

`price_shifted_antipodal_capacity.py` retains both the original pair
charge and all full/exterior/lattice count-budget charges, then changes
only the directed exterior capacity terms. The old frozen dual is
`nested_conditional_jointcount_fullgauge_dual_geometry.npz`, SHA-256
`9608d101717d76bb61a4b2ece9bb80c61598c69886cbc354dab9cd0266c9ac27`.

The ordinary antipodal energy remains exactly zero on lattice rows by
Section 2; the reduction of the former count energy is still charged.

| Profile | Shift delta | Raw finite-dual price |
|---|---:|---:|
| Inner quartic B4 | -0.0043875 | -0.9805052 |
| Inner quartic B0 | 0.1954678 | -0.05916385 |
| Indicator through R26 | 0.1955367 | -0.05700112 |
| Indicator through R25 | 0.1955374 | -0.05278120 |

At zero shift the B0 price was positive, about `1.44774`. Thus allowing
the test center to move produces additional negative directions without
losing central-symmetry equality. The search uses a finite grid and a
bounded numerical refinement; it does not prove an optimal shift.
These numbers are raw reduced costs, with raw objective equal to eight
times density. No nonzero feasible allocation coefficient has been
established, and the table is not an upper bound on packing density.

This note and the pricing formulas were submitted for independent
analytic audit before use in any global certificate. No claim of
novelty for the movable-center construction is made.
