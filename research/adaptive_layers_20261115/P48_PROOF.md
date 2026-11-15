# A certified deep binary plane above P48q

Three nonzero half-cosets form a binary plane. A fresh complete enumeration
excludes all numerator vectors of squared norm <=12 in each of these
classes, so their half-coset squared depths are at least 4. Extend the
plane by an anisotropic two-space: all twelve additional classes are
nonsingular and have numerator norm >10, hence half-coset depth >=7/2.
The binary form in this basis is `q(a)=a2+a3+a2*a3`, with zero-based indices.

All fifteen positive enumerations completed, totaling **1,002,136,588
nodes**, and an exact-witness negative control passed. Inputs, completed
reports, compiler flags and hashes are archived in [p48/replay](p48/replay/).
This computation uses the documented outward-rounded IEEE-754 arithmetic
assumptions of the existing enumerator. P48q's published minimum 6 is inherited.

With the first r rows B of the explicit binary basis, use vertical Grams
from the leading blocks of

```
M4 = [[ 2, -1,    0,    0],
      [-1,  2,   -1,   -1],
      [ 0, -1,  5/2, -1/4],
      [ 0, -1, -1/4,  5/2]].
```

The full Gram is

```
[[G, G*B^T/2],
 [B*G/2, B*G*B^T/4 + M]].
```

Its determinant is det(M). All-even nonzero vertical coefficients have
squared norm at least 6; the three deep parities have norm at least 2;
all nonsingular parities have norm at least 5/2. Adding the horizontal
bounds gives minimum 6. Positivity, determinants and every potentially
violating vertical vector are checked exactly using inverse-metric
coefficient bounds.

| Dimension | det(M) | Center-density square | Gain over preceding bound |
| --- | ---: | --- | ---: |
| 50 | 3 | `(3/2)^50/3` | 15.470054% |
| 51 | 11/2 | `(3/2)^51/(11/2)` | 11.803399% |
| 52 | 121/16 | `(3/2)^52/(121/16)` | 11.340443% |

The explicit residues and full Gram formula define the packings. Horizontal
representatives are supplied in [p48/certificate.npz](p48/certificate.npz).

[Exact setup and saved-report verifier](verify_p48.py) · [Results](p48_results.json).
A fresh complete replay is available with `--replay NEW_DIRECTORY`.
