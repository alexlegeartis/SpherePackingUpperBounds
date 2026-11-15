# Exact recycling of pair slack by radial count gauges

This note gives a locally nonnegative kernel construction and its finite-dual
pricing rule. It does not establish a new numerical sphere-packing bound.
The accepted certificate at `pointed/nested_repair2` is unchanged.

## 1. Local convention and the recycling theorem

Let `C` be a finite set of nonzero neighbors of an anchor at the origin. Assume
that `{0} union C` has all distances at least one, and retain only neighbors
whose radii are in `[1,R]`. A local kernel is written as

\[
E_0(C)=\kappa+2\sum_{x\in C}h(|x|)
       +\sum_{x,y\in C}K(x,y)\ge0.
\]

The double sum is ordered and includes its diagonal. The scalar auxiliary
function is `f`, and the complete base pair expression is

\[
q(r)=f(r)+2h(r)+K(x,x),\qquad |x|=r.
\]

Here `h` and `K` include every base component, including integer count
quadratics, pointed probes, tents, and PSD kernels. Replacing this expression
by `f` plus the PSD diagonal is justified only on a region where all other
base pair contributions vanish.

Partition `[1,R]` into finitely many radial regions. For each pivot region
`i`, choose a nonnegative radial gauge `psi_i` such that

\[
\sum_{x\in C}\psi_i(|x|)\le1
\quad\hbox{for every admissible neighborhood }C,
\qquad 0\le\psi_i(r)\le1-\epsilon
\]

for some fixed positive `epsilon`. Gauges and weights are zero outside the
common support. Different pivot regions may use different valid gauges.
Assume `q(r)<=0` throughout the support and define

\[
w(r)=\frac{-q(r)}{1-\psi_{i(r)}(r)}\ge0.
\]

Then the local energy

\[
E_G(C)=\sum_{y\in C}w(|y|)
 \left[1-\sum_{x\in C}\psi_{i(|y|)}(|x|)\right]
\tag{1}
\]

is nonnegative: every summand is nonnegative separately. Its symmetric kernel
encoding is

\[
\kappa_G=0,\qquad h_G(r)=\frac12w(r),
\]

\[
K_G(x,y)=-\frac12\left[
 w(|y|)\psi_{i(|y|)}(|x|)
 +w(|x|)\psi_{i(|x|)}(|y|)\right].
\tag{2}
\]

Summing (2) over ordered pairs gives exactly (1). In particular,

\[
K_G(x,x)=-w(r)\psi_{i(r)}(r),\qquad
2h_G(r)+K_G(x,x)=w(r)[1-\psi_{i(r)}(r)]=-q(r).
\tag{3}
\]

Thus the new actual pair expression is identically zero. There is no new
origin cost: the entire old `kappa`, including constants from count
polynomials, is retained without alteration. The construction uses local
nonnegativity rather than an unstated PSD condition.

Every nonorigin entry in (2) is nonpositive. Consequently every full cyclic
triangle expression, and every exterior expression with only one surviving
kernel pair, can only decrease. The Fourier condition and scalar tail of
`f` are unchanged. In the established local-energy upper-bound theorem, this
transforms any feasible base certificate into another feasible certificate
with the same objective and pair equality throughout `[1,R]`.

No assertion here supplies a feasible base certificate at a desired
objective. During numerical discovery the inequalities `q<=0` remain
essential: clipping a positive numerical residual to zero would destroy the
affine row model.

All endpoint values use the actual closed radial counts. Distinct limits at
a jump can also be checked as separate closure branches. Weights and both
nonzero arguments of the kernel are truncated to `[1,R]`; the anchor is never
included in `C`. These conventions also cover empty neighborhoods, for which
the added local energy and origin constant are both zero.

## 2. Three exact gauges from the accepted count implications

In increasing order, the current five radii are

\[
(r_0,r_1,r_2,r_3,r_4)=
\left(\frac{501473}{500000},\frac{1005643}{1000000},
\frac{50379}{50000},\frac{50767}{50000},\frac{102597}{100000}\right).
\]

Write `N_j` for the number of neighbors of radius at most `r_j`. The exact
count envelope gives `N_2<=24`, `N_3<=25`, and `N_4<=26`. The exact affine
latitude certificates give

\[
N_0=24\ \Longrightarrow\ N_4=24,
\qquad
N_1=24\ \Longrightarrow\ N_3=24.
\tag{4}
\]

Their proof data are in `spherical_sdp/accepted_count_envelope.json` and
`conditional/affine_latitude_certificates.json`. This note takes those
verified geometric implications as premises.

The following three gauges therefore satisfy the capacity inequality:

\[
\psi_{24}=\frac{\mathbf1_{[1,r_2]}}{24},\qquad
\psi_{25}=\frac{\mathbf1_{[1,r_1]}+\mathbf1_{[1,r_3]}}{48},\qquad
\psi_{26}=\frac{2\mathbf1_{[1,r_0]}+\mathbf1_{[1,r_4]}}{72}.
\tag{5}
\]

For example, if `N_1=24`, (4) bounds `N_1+N_3` by 48; otherwise it is at most
`23+25=48`. Similarly, `2N_0+N_4<=72` follows from the two alternatives
`N_0=24` and `N_0<=23`. Every gauge in (5) has value `1/24` at contact and is
at most `1/24` everywhere. Hence all denominators in (3) are at least `23/24`.

### Pareto reduction

Consider gauges

\[
\psi=\sum_{j=0}^4a_j\mathbf1_{[1,r_j]},\qquad
a_j\ge0,\quad \sum_j a_j=\frac1{24},
\]

with

\[
a_0\ge2a_4,\qquad a_0+a_1\ge a_3+2a_4.
\tag{6}
\]

These conditions imply capacity at most one. Indeed, the three count cases
`N_0=24`, `N_0<=23,N_1=24`, and `N_1<=23` bound the weighted counts by,
respectively,

\[
1,\qquad 1-a_0+2a_4,\qquad
1-a_0-a_1+a_3+2a_4.
\]

Every gauge satisfying (6) is pointwise dominated by the gauge with

\[
(a'_0,a'_1,a'_2,a'_3,a'_4)=
\left(2a_4,a_3,\frac1{24}-2a_3-3a_4,a_3,a_4\right).
\tag{7}
\]

Nonnegativity of its middle coefficient follows from (6). The increases on
the two affected bands are `a_0-2a_4` and
`a_0+a_1-a_3-2a_4`; all other band values agree. Expression (7) is the convex
combination of (5) with coefficients

\[
\theta_{24}=1-48a_3-72a_4,\qquad
\theta_{25}=48a_3,\qquad\theta_{26}=72a_4.
\]

Increasing a gauge pointwise makes every directed negative kernel term in
(2) at least as strong for fixed `q<=0`. It also decreases the local energy
at every neighborhood: for a pivot `y`, the ratio
`(1-sum_C psi)/(1-psi(y))` decreases when any gauge value increases. Thus the
three-vertex hull contains a pointwise dominating representative for every
gauge in (6).

## 3. Cyclic coefficients and finite-dual pricing

Set `g(r)=-q(r)>=0`. A directed target--pivot occurrence `(x,y)` in a cyclic
row contributes

\[
-\frac12\frac{g(|y|)\psi_{i(|y|)}(|x|)}
                    {1-\psi_{i(|y|)}(|y|)}.
\tag{8}
\]

A full triangle has six such directed occurrences, two for each of its
three kernel pairs. An exterior row has only two directed occurrences.
For the basic gauge `psi_24`, a full triangle with exactly one side in the
contact shell receives `(q(b)+q(c))/48`, where `b,c` are its other sides.
With two shell sides `a,b`, it receives `(q(a)+q(b))/46+q(c)/24`; with three,
it receives `(q(a)+q(b)+q(c))/23`. These are changes in the cyclic kernel
expression, not changes in density.

For a finite-dual diagnostic, multiply (8) by the corresponding raw cyclic
row multiplier and sum. A lattice guard with shell radii `rho_m`, shell
counts `n_m`, and multiplier `eta` adds, for a pivot shell `n` in region `i`,

\[
\eta n_n g(\rho_n)
\frac{1-\sum_m n_m\psi_i(\rho_m)}{1-\psi_i(\rho_n)}.
\tag{9}
\]

Equations (8)--(9) include the complete new local energy and have no hidden
origin constant. Existing nonnegative-pair-slack rows are retained to prove
`q<=0`; they are not reinterpreted as actual nonzero pair terms after the
lift. Replacing an older slack budget requires subtracting that older
budget's cyclic and lattice changes first.

### Why three choices per pivot region suffice for pricing

Assume each candidate gauge is constant on every fixed radial region. For
a pivot region `i`, write `psi_i=sum_j theta_j psi_j`, with nonnegative
`theta_j` summing to one. For fixed base `q` and fixed dual multipliers, its
total contribution from (8)--(9) is

\[
F_i(\theta)=\frac{\sum_j\theta_j A_{ij}}
                  {\sum_j\theta_j D_{ij}},
\qquad D_{ij}=1-\psi_j(\hbox{region }i)>0.
\]

This is a denominator-weighted average of the vertex ratios:

\[
F_i(\theta)=\sum_j
\frac{\theta_jD_{ij}}{\sum_k\theta_kD_{ik}}
\frac{A_{ij}}{D_{ij}}.
\]

Its minimum is attained at a vertex. Contributions from different pivot
regions are independent, so a minimum of this fixed-dual pricing problem
can be obtained by choosing the best of the three gauges separately for
each region. If a candidate gauge has a jump inside a region, the partition
must be refined before using this argument.

This pricing conclusion is narrower than a claim that pure vertex choices
always give the best globally feasible certificate. A mixture may help
meet several competing constraints. For fixed `q`, all mixture kernels in
a region are exactly the convex hull of its three vertex lift kernels,
using the denominator-weighted coefficients above. Joint optimization of
the base coefficients and gauge mixtures is not made linear merely by this
vertex pricing calculation.

The current raw master objective is eight times the variable part of
density. Therefore raw reduced costs and dual-weighted row changes must be
divided by eight to express density units. Even after that conversion, a
finite-dual row effect is not an attained objective change or a packing
bound.

## 4. Redundancy and allocations of slack

Fix a gauge `psi`. An arbitrary nonnegative weight `u` gives the valid local
product `(sum_C u)(1-sum_C psi)`. Before the full lift it adds

\[
\delta h=u/2,\quad
\delta K=-\operatorname{sym}(u\psi)/2,\quad
\delta q=(1-\psi)u.
\]

Applying the same full lift afterwards decreases its residual weight by
exactly `u`. Thus these changes cancel identically in the final kernel.
The product ray can tighten the retained base condition `q<=0`, but cannot
enlarge the final-kernel family; setting its coefficient to zero is always
at least as feasible. This includes squared weights `u=p^2` and explains why
same-gauge count products should not be added as new variables after the
lift.

Products using a different gauge can remain useful. A broader exact
construction splits the available slack into nonnegative parts

\[
g(r)=\sum_j g_j(r),\qquad g_j(r)\ge0,
\]

and uses weight `w_j=g_j/(1-psi_j)` with gauge `psi_j`. The added pair costs
sum to `g`; every local product remains nonnegative. If the `g_j` are
represented directly, this is a convex allocation of pair slack, subject
to their pointwise nonnegativity and the exact sum constraint. Choosing
fixed fractions by pivot region is a small special case. A usable numerical
implementation must still certify those radial nonnegativity conditions.

## 5. Independent scope check

The independent certificate-audit agent checked the local energy,
symmetrization, diagonal, zero origin cost, and vertex argument. Its one
required scope condition was explicit constancy of every candidate gauge
on the pivot region; that condition is included above. All geometric input
is confined to the already verified count and conditional-count premises.
No assertion uses a property specific to D4 to infer validity for arbitrary
packings.
