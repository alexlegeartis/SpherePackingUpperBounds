# Broad antipodal graph: capacity 1208 and the saving 3/8000000

Family: all sphere packings in R^4. This note proves a bound on one graph inside an anchored neighbourhood. With the accounting of the base certificate it lowers the density bound by exactly 3/8000000.

## Setting

Place a packing centre at the origin. Its neighbours are the other centres y with |y| <= 7/4. All centres are at mutual distance at least 1.

The broad antipodal graph has these vertices and edges.

- Vertices: neighbours with radius in J = [33/28, 7/4].
- Edges: pairs of distinct vertices y, z whose directional cosine t = (y . z)/(|y| |z|) lies in [-1, -7/8].

Both intervals are closed. The edge count e counts ordered pairs, so each unordered edge counts twice.

The base certificate (stage 1, `data/v7/graph_energies.json`) carries the energy delta (1358 - e) with delta = 1/50000000. Its constant 1358 delta is paid in the origin constant, so it sits inside U_raw. The census floors of stages 2 and 3 use only the fact that this energy is nonnegative.

## Theorem

Every anchored neighbourhood has e <= 1208. So delta (1358 - e) >= 150 delta = 3/1000000. Since the scalar integral is pi^2/4 and the density is (f(0) + C)/8 in the convention of the base certificate, this extra energy lowers the density bound by 3/1000000 / 8 = 3/8000000.

## Step 1. Directions in the annulus are separated

Write y = r u and z = s v with unit vectors u, v. Unit separation gives

    u . v <= h(r, s) = (r^2 + s^2 - 1)/(2 r s).

For fixed s the derivative in r has numerator r^2 - s^2 + 1, which increases in r. Any interior critical point is a minimum, so the maximum over the square J x J is at a corner. With L = 33/28 and R = 7/4 the corner values are

    h(L, L) = 697/1089,   h(L, R) = h(R, L) = h(R, R) = 41/49.

So all distinct directions satisfy u . v <= rho = 41/49. The lower end 33/28 equals R - 1/R. That choice keeps the mixed corner at the same value as the outer corner.

## Step 2. At most 151 vertices (exact certificate)

On S^3 use the normalized zonal polynomials Z_l(t) = U_l(t)/(l + 1), where U_l is the Chebyshev polynomial of the second kind. The Gegenbauer parameter in dimension 4 is 1. For any finite list of unit vectors the addition formula gives sum over i, j of Z_l(u_i . u_j) >= 0.

The certificate `src/annular_count_151_certificate.json` gives a degree-20 polynomial p = sum a_l Z_l with rational coefficients. It has a_0 = 999999/1000000 > 0 and a_l >= 0 for l >= 1. The checker proves p(t) <= 0 on the whole closed interval [-1, 41/49]. It rebuilds the monomial form of p, converts it to Bernstein form and subdivides with exact rational de Casteljau steps. The proof uses 93 nodes, 47 leaves and depth 12.

For n vertices, the addition formula and the off-diagonal sign give

    a_0 n^2 <= sum over i, j of p(u_i . u_j) <= n p(1).

So n <= p(1)/a_0 = 151196972764393/999999000000 < 152, and n <= 151.

## Step 3. Projection at a vertex

Fix a vertex with direction u and put w = -u. Every neighbour v of u in the graph has h = v . w >= beta = 7/8.

If some neighbour has h = 1, then v = w. A second neighbour v' would have v . v' >= 7/8 > 41/49, which contradicts Step 1. So that vertex has degree at most 1.

Otherwise write v_i = h_i w + sqrt(1 - h_i^2) x_i with unit vectors x_i perpendicular to w. These x_i are unit vectors in a copy of R^3. From v_i . v_j <= rho,

    x_i . x_j <= Q(h_i, h_j) = (rho - h_i h_j)/sqrt((1 - h_i^2)(1 - h_j^2)).

The derivative of Q in h_i has the sign of rho h_i - h_j, which is at most rho - beta < 0 on [beta, 1). So Q is largest at h_i = h_j = beta, and

    x_i . x_j <= Q(7/8, 7/8) = (41/49 - 49/64)/(1 - 49/64) = 223/735.

## Step 4. Nine unit vectors in R^3 (exact certificate)

The shipped certificate `data/projected_eight.json` is a nonnegative
combination of normalized Legendre polynomials. The checker
`src/verify_projected_eight.py` rebuilds those polynomials by their exact
recurrence, maps `[-1,223/735]` to `[0,1]`, and proves the sign of the
negative polynomial with rational Bernstein subdivision. It checks the
endpoint separately. The exact objective is

    89804130521/9999900000 < 9.

The addition formula therefore gives at most eight unit vectors in R^3 with
all pairwise inner products at most `223/735`. The classical Schütte--van der
Waerden nine-point theorem remains a cited alternative, but is not required by
this replay.

## Step 5. Degree at most 8, edges at most 1208

Nine neighbours of one vertex would give nine unit vectors x_i in R^3 with all inner products at most 223/735. The exact certificate in Step 4 forbids this. The geometric gap `1/3 - 223/735 = 22/735 > 0` is also checked by `gap_check.py`. Every vertex has degree at most 8 (or at most 1 in the case of Step 3 with h = 1). The ordered edge count is the degree sum, so

    e <= 151 * 8 = 1208.

The older bound 1358 used a degree-nine polynomial bound and parity. This argument does not use that polynomial.

## Step 6. The edges are exterior

A selected pair has |y - z|^2 = r^2 + s^2 - 2 r s t >= (2 + 7/4) L^2 = 16335/3136 > 49/16. So |y - z| > 7/4. The graph energy has no diagonal term and no coupling to the origin. It enters only through exterior triangles, as in the base certificate.

## Step 7. The density saving

The graph charges of the base certificate, as (weight, capacity), are (1/50000000, 1358), (1/5000000, 120), (1/5000000, 604) and (1/100000000, 1880). Their sum is 4769/25000000. The base accounting is

    U_raw = objective(origin.json) + 4769/25000000 / 8.

`run.py` checks this identity exactly against `census_floors.json` and the stage 2 acceptance file.

The census floors bound the anchored energy E from below using only delta (1358 - e) >= 0 for this graph. By the theorem, that term is at least 150 delta = 3/1000000 at every anchor. So every floor rises by 3/1000000, and the bound

    density <= U_raw - a/8

becomes U_raw - a/8 - 3/8000000. Equivalently, the capacity 1358 in U_raw may be replaced by 1208. This is the term `graph_saving` = 3/8000000 in `census_floors.json` and in stage 5.

## Evidence

- Steps 1, 3, 5, 6 and 7: written proof plus exact rational checks in `run.py` and `src/gap_check.py`.
- Step 2: exact certificate, rechecked by `src/verify_annular_count_151.py` (no floating point).
- Step 4: exact rational Legendre/Bernstein replay from `data/projected_eight.json`; the classical theorem is attribution only.
