# Current56D input and elementary propagated57D/58D benchmarks

## Published input

The current Henry Cohn table gives Schuett's56-dimensional centre density

```
d = 5^28/(2^28 * 29^4) = 196213.21381197087...
```

Source: <https://cohn.mit.edu/sphere-packing/>, row56, referring to
M. Schuett and T. Shioda, *Mordell-Weil Lattices*, Springer2019, p400.
The primary catalogue independently lists `2^17.58` and the same reference:
<https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html>.
Its page header still says last modified February2012, but its actual56D
row includes the2019 result; the header is not a reliable date for every row.

The following bounds are mathematical consequences of that published input.
No claim is made that their dimension-specific numerical values were
separately published, nor that this is an exhaustive envelope of all known
constructions.

## Greedy torus offsets

Normalize a parent lattice L in dimension n to minimum4. If its centre
density is d, then its covolume is1/d. Let V_n be the volume of the unit
ball. For any s>=2, one can choose s points h_1,...,h_s in the torus R^n/L
with pairwise torus distances at least R whenever

```
(s-1) d V_n R^n < 1.
```

Indeed, after selecting j<s points, the union of their forbidden open
radius-R balls has torus volume at most j V_n R^n, strictly less than1/d.
Choose the next point outside it. The image of a ball in the torus cannot
have greater volume than its Euclidean volume, so no injectivity assumption
is needed. Compactness yields the non-strict limiting separation as well.

## A general colored-layer construction

Let K be a rank-k lattice normalized to minimum1. Suppose it has a periodic
s-coloring in which distinct points of the same color have squared
separation at least lambda. At every vertical point y in tK place a full
horizontal translate L+h_color(y). Choose

```
t^2 = max(4-R^2, 4/lambda).
```

Horizontal points within the same layer are separated by at least2.
Points in different colors have squared separation at least R^2+t^2>=4.
Points of the same color in distinct layers have vertical squared
separation at least lambda*t^2>=4. This gives a periodic packing of
unit-radius balls, with centre density

```
d_new = d / (t^k covol(K))
      = d * 2^k * delta(K) / (t^2)^(k/2).
```

The torus offsets need not form a subgroup, so the result need not be a
lattice. Each layer has exactly the same density d.

For a general easy family, color K by K/2K: then s=2^k and lambda=4.
For k=2, a better choice is the triangular lattice with its proper
three-coloring: s=3, lambda=3, and covol(K)=sqrt(3)/2.

## Consequences in57D and58D

For57D take K=Z, two colors, lambda=4. The limiting formula is

```
R^2 = (1/(d V_56))^(1/28),
delta_57 >= d/sqrt(4-R^2) = 151706.1088095999...
```

For58D take the three-colored triangular lattice. Its same-color
sublattice has minimum3. The limiting formula is

```
R^2 = (1/(2d V_56))^(1/28),
delta_58 >= 2d/(sqrt(3)*(4-R^2)) = 130984.5416742293...
```

The same-color condition is satisfied in each case. To obtain completely
algebraic bounds with strict volume inequalities, the checker instead
uses these rational squared separations:

|Dimension|R^2|t^2|Certified centre density|
|---:|---:|---:|---:|
|57|2327/1000|1673/1000|d sqrt(1000/1673) =151698.1486935438...|
|58|227/100|173/100|200d/(173 sqrt(3)) =130963.8749282115...|

`verify_torus_layer_bounds.py` verifies the required union-volume
inequalities as exact rational inequalities using pi<355/113. It also
certifies that pi bound from Machin's identity and alternating arctangent
series. Its output records exact squared centre densities.

These comparisons already exceed the Ne64 bare-section candidates of
57530.1758... in57D and115060.3516... in58D. The later exact Gamma72/Leech
section constructed by the other agent gives an even stronger58D
comparison, namely2^18/sqrt(3)=151348.9089....
