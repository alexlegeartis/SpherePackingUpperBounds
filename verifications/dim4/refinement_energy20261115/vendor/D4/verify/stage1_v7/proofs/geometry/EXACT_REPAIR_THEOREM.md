# Three exact graph energies for the recovered D4 candidate

Status: new exact **local geometric theorems**. No uniform bound on a base
candidate's cyclic expression is proved by these files. No global density
bound is claimed here.

Normalize packing-center separation to one. Fix a center at the origin.
Every interval in the following definitions is closed except the indicated
open endpoint. All graph edge counts are ordered and hence even.

## 1. Three count bounds with no inherited prerequisites

| Graph | Vertex classes | Selected edges | Exact ordered-edge bound |
|---|---|---|---:|
| A | First class `[1,103/100]`; second class `(103/100,3/2]` | Cross-class directional inner product in `[-3/4,-1/2]` | 820 |
| B | First class `[1,103/100]`; second class `[3/2,7/4]` | Cross-class directional inner product in `[-19/20,-3/4]` | 1102 |
| C | One class `[4/3,7/4]` | Mutual distance in `[4/3,7/4]` | 9630 |

These theorems need no count24, count26, count96, count155, count227, or
previous global packing certificate. The proofs use only separation,
the spherical harmonic addition formula, explicit rational positive
Gram rays, and elementary exact polynomial and integer calculations.

For A and B, use the normalized four-dimensional zonals
`Z_l(t)=U_l(t)/(l+1)` and a symmetric matrix polynomial
`F(t)=sum_l A_l Z_l(t)`. Every nonconstant `A_l` is a stated nonnegative
rational sum of matrices `v v^T`. The diagonal polynomials are nonpositive
on the entire possible same-class angular domain. The cross polynomial
is nonpositive on its entire domain and at most minus one on the selected
angular band. All these complete interval claims are proved by rational
Bernstein subdivision.

For any radial intervals, the maximum of
`(r^2+s^2-1)/(2rs)` over the radial rectangle occurs at a corner: its
derivative in either variable has a numerator that crosses zero only from
negative to positive, so an interior critical point is a minimum. This
gives the exact angular domains used by the checker, including closures.

Writing class counts as `n,m`, harmonic positivity and the graph sign
conditions give

`e <= n F_NN(1) + m F_MM(1) - (n,m) A_0 (n,m)^T`.

The two constant matrices `A_0` are positive definite. The checker
completes their quadratic squares and proves that every integer `n>=100`
gives a negative upper bound even after maximizing over all real `m`.
For each `0<=n<100`, it checks zero and both nearest integers to the
quadratic vertex in `m`. Thus the maximization covers **all** nonnegative
integer counts. The exact maxima are approximately `821.7156123828156`
and `1103.6483818488912`; parity gives 820 and 1102.

For C, a degree-26 polynomial with nonnegative normalized-U coefficients
is nonpositive on `[-1,41/49]`, and at most minus one on
`[71/512,313/441]`. These are exactly the angular bounds needed for
radii and selected distances in `[4/3,7/4]`. The concave bound
`e <= n f(1)-f_0 n^2` has integer maximum approximately
`9631.774453230664`, giving 9630 by parity. Again, no independent
upper bound on `n` is needed.

## 2. Exact local energies and where they act

For arbitrary nonnegative rational `delta_A,delta_B,delta_C`, each
`E_j=delta_j(C_j-e_j)` is a nonnegative local energy. Its origin constant
is `delta_j C_j`; its linear term is zero; and its nonorigin kernel is
minus `delta_j` times its symmetric graph-edge indicator.

Every diagonal value is zero. For A and B, a self-pair has inner product
one and cannot be selected. For C, a self-pair has mutual distance zero.
Consequently these corrections leave every pair inequality and any
existing complete pair-count gauge unchanged. The new constants are
separate local constants, not modifications of a retained PSD completion.

For A, every selected mutual distance is strictly greater than `7/4`:
its square is at least `1+(103/100)^2+103/100=30909/10000>49/16`.
For B, the corresponding square is at least `1+9/4+9/4=11/2>49/16`.
Thus A and B affect only exterior triangles when active support is `7/4`.
They subtract one `delta_j` at the active vertex and zero at the two
vertices whose incident edge exceeds support.

For C, a full triangle with every side in `[4/3,7/4]` is selected at
each of its three vertices. Its cyclic decrement is exactly `3 delta_C`.
Every other change from any of the three corrections is nonpositive.

## 3. A concrete conditional budget

Choose

`delta_A=1/1000000`, `delta_B=3/10000000`, `delta_C=1/15000000`.

The total density cost in the retained normalization is exactly

`(820 delta_A + 1102 delta_B + 9630 delta_C)/8 = 8963/40000000 = .000224075`.

Added to the recovered candidate's recorded strict threshold
`629973666871/10^12`, this yields the **conditional** threshold

`630197741871/10^12 = .630197741871`.

The unchanged base candidate would need its complete cyclic expression
bounded by `delta_A` on graph-A exterior support, by `delta_B` on
graph-B exterior support, by `2e-7` on full triangles with all sides
in `[4/3,7/4]`, and by zero elsewhere. Where selected graph supports
overlap, the allowed corrections add. Every other theorem hypothesis
must hold as well.

These weights repair every positive witness family in the recovered
`alpha0015_new_separator_exact_points` record. Pointwise witness repair
does **not** imply a uniform inequality on any continuous branch.

`graph_repair_bounds.py` supplies an exact conservative lower bound on
the cyclic decrement over a closed interval box, with the same exterior
coordinate as the retained verifier. A box crossing a support boundary
gets zero from that graph until further subdivision. It has no global
acceptance path.

## 4. Reproduction

Run `python verify_exterior_graphs.py` and
`python verify_long_triangle_graph.py` in this directory. Python's
standard library is sufficient. The numerical discovery/rationalization
scripts are not needed for either exact replay.
