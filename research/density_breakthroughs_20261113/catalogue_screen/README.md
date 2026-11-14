# The 56D contact-preserving pencil: exact exclusion

**No new density record.** The entire one-parameter family discovered here has
center density at most `(3/2)^28 = 85222.69299239293`, below the reference
`(5/2)^28 / 29^4 = 196213.21381197087`.

The source is the Nebe–Sloane catalogue entry
[L_56a](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/L_56a.html),
archived in `sources/`. Its integral Gram matrix G has determinant `29^28`
and stated minimum 28. Exact group orbits produced 53,940 antipodal norm 28
vectors, matching the catalogue's 107,880 signed contacts. Their quadratic
constraints have numerical rank 1595 out of 1596. An integer symmetric matrix D
was reconstructed and checked to vanish on every one of those contacts.

The exact identities are

```
tr(G^-1 D) = 0,
(G^-1 D)^2 = I/29,
det(G+tD) = (29-t^2)^28.
```

The positive-definite domain is `|t| < sqrt(29)`. Although determinant descent
is available, it does not give a new packing record.

## A global argument within the pencil

The compact file `L56a_pencil_certificate.npz` contains G,D, two integral
coordinate changes U,R, and two integer vectors v,w. The exact verifier checks
positive definiteness of G by rational LDL, determinant and integrality, and

```
det(U)=1, det(R)=+/-1,
U^T G U = (27 G+145 D)/2,
U^T D U = (5 G+27 D)/2,
R^T G R = G, R^T D R = -D,
v^T G v=28, v^T D v=0,
w^T G w=32, w^T D w=-4.
```

Thus reflection changes t to -t. The inverse of U, followed by positive
rescaling, changes t to

```
g(t)=(27t-145)/(27-5t).
```

For `5<t<sqrt(29)`, the denominator is positive and `-5<g(t)<t`. If an
iteration remained above 5 forever, its decreasing limit would be a fixed
point of g. The only fixed points are `+/-sqrt(29)`, impossible below the
initial parameter. Every member of the positive-definite pencil is therefore
equivalent, up to unimodular change of basis and scale, to one with `0<=t<=5`.

The two displayed vectors imply that the true minimum mu(t) satisfies
`mu(t)<=28` and `mu(t)<=32-4t`. Consequently the center density is at most

```
[49/(29-t^2)]^14                      for 0<=t<=1,
[(8-t)^2/(29-t^2)]^14                 for 1<=t<=5.
```

The first expression increases to `(7/4)^14`. In the second, the derivative
changes from negative to positive at `t=29/8`, so its maximum is at an
endpoint. The larger endpoint value is `(9/4)^14=(3/2)^28` at t=5. This bound
does not require a complete shortest-vector enumeration, and it permits
changes in the minimum away from the original contact shell.

Run from the repository root:

```bash
python research/density_breakthroughs_20261113/catalogue_screen/verify56_pencil.py
```

This needs only NumPy and SymPy and the compact certificate. `orbit56.py`,
`contact56.py`, `flex56.py`, and `test56.py` record the discovery route;
the large orbit and shell caches are reproducible and excluded from git.
The other archived catalogue pages are a metadata screen, not evidence of
successful deformations in those dimensions.
