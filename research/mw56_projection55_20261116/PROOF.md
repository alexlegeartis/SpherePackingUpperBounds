# MW56 projections against the full Leterrier55 target

The acceptance target in 55D is **116,783.149646304716…**, the full
center density printed in Leterrier's 2024 thesis, Proposition 3.2.22.
Our existing certified 91,333.867143595418… is below it. Improving a
weaker repaired construction does not clear this target.

This search starts from the published **56D Mordell–Weil lattice** of
minimum 10 and Gram determinant 29^8, then projects and combines cosets
in 55D. Five proposals had densities above the target *if* their minimum
were 10. All five have exact shorter-vector counterexamples. A stronger
certificate excludes **every subset of projected cosets, at every positive
radius, in each of the 24 supplied projection directions**.

**No new packing lower bound or SOTA is claimed.** The largest certified
density ceiling among these families is **103,140.241156608598…**,
below the target. This is not a global 55D upper bound and does not exclude
other MW56 projection directions, arbitrary shifts, larger periods, or
other constructions.

## External targets and starting lattice

[targets.json](targets.json) records the full named targets for 54D–56D:

| Dimension | External work | Squared center-density target |
| --- | --- | --- |
| 54 | Leterrier (2022), Corollary 1.4 | `(5/2)^54 / 3^25` |
| 55 | Leterrier (2024), Proposition 3.2.22, printed claim | `(5/2)^55 / (2*3^24)` |
| 56 | Shioda's MW56; Schütt–Shioda, Theorem 13.35 at p=29 | `(5/2)^56 / 29^8` |

The original [2022 paper](https://doi.org/10.1007/s40993-022-00321-0)
establishes the 54D lattice. The
[2024 thesis](https://infoscience.epfl.ch/record/307330/files/EPFL_TH10417.pdf)
prints the larger 55D claim, but its displayed lattice contains a norm-6
period vector; see the [exact correction](../tilted55_20261116/PROOF.md).
That defect does not prove the numerical density unattainable. We retain
the full number as the challenge target and do not use it as a proved
lower bound in product or section propagation.

The valid strong parent for this experiment is MW56 from
[Schütt–Shioda, Theorem 13.35, p.357](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf).
The verifier reconstructs its Gram matrix G using the published polynomial
sections and height pairing through [the existing source audit](../cm58_20261116/verify.py).
It checks `det G = 29^8` and that `H = 29 G^-1` is integral. Source locations
and the dated review are recorded in [sources.json](sources.json).

## Projected cosets

All coordinates below are row vectors relative to the parent basis.
For each supplied primitive integer row b, let the dual normal have
coordinates `alpha = b G^-1`. The certificate checks

```
||alpha||² = m/29,
v0 = 29 alpha = b H is primitive in the parent lattice,
v0 b^T = m,
||v0||² = 29m.
```

Let S be the rank-55 section defined by `v b^T = 0`, and let pi be
orthogonal projection onto the normal's perpendicular space. Exact Hermite
normal form constructs a basis of S and verifies

```
det Gram(S) = 29^7 m.
```

Primitivity of b makes its integer label map onto Z. Primitivity of v0
makes the kernel of projection on the parent lattice exactly Z v0.
Consequently `pi(L)/S` is cyclic of order m. A projected coset has label
`j = v b^T mod m`. Its vectors have squared norm

```
||pi(v)||² = v G v^T - 29 (v b^T)^2 / m.
```

The 24 supplied directions have m between 98 and 128. They came from
reduction of the dual Gram form; they are not a complete classification
of short dual normals or of their automorphism orbits.

## Proposals large enough to matter, and exact obstructions

For M distinct cosets of S, a squared minimum of 10 would give

```
delta² = M² (5/2)^55 / (29^7 m).
```

Here are five explicit proposals with this nominal density above the full
printed target. The case numbers are zero-based indices in the certificate.
Coset labels are elements of Z/mZ.

| Case | m | Chosen labels | Density assuming minimum 10 | Exact shorter squared norm |
| --- | --- | --- | ---: | ---: |
| 1 | 98 | 0, 28 | 135,012.569813 | 8 |
| 5 | 104 | 0, 52 | 131,060.123194 | 8 |
| 8 | 108 | 0, 36, 72 | 192,915.277860 | 8 |
| 14 | 116 | 0, 4 | 124,096.132533 | 8 |
| 17 | 120 | 0, 60 | 122,010.334697 | 6 |

These are rejected proposals, **not achieved densities**. The witness
at the indicated nonzero label belongs to a difference of two chosen
cosets. The verifier evaluates every witness by exact integer arithmetic
and checks that even the density ceiling obtained by rescaling to its
shorter distance is below the target.

Discovery used lattice reduction and bounded short-vector walks, followed
by targeted affine-coset enumeration and reduction of index-two and
index-three overlattices. Completeness of those searches is not needed:
the proof uses the explicit resulting vectors, not a search timeout or
failure to find a vector.

## Excluding every coset subset and every radius in these families

[certificate.npz](certificate.npz) contains **2,704 parent vectors**,
one for each residue of each family. Their projected squared norms give
upper bounds U_j on the corresponding distances. The exact checks include
`U_j = U_(-j) > 0` and `U_0 = 10`. The residue-zero witness is nonzero
after projection, so every nonempty union has minimum at most 10.

For any proposed positive squared minimum rho, a necessary condition on
two chosen labels is `U_(i-j) >= rho`. Form the undirected graph on Z/mZ
with these allowed differences. Every admissible subset is a clique.
Translation lets us assume a nonempty clique contains 0. The other
vertices are then in the allowed-neighbor graph of 0. If that induced
graph has a proper coloring with c colors, its cliques have at most c
vertices, so the full subset has at most `M = c+1` translates.

[colorings.json](colorings.json) supplies color classes as integer bitmasks.
The verifier checks that they partition the allowed neighbors and contain
no allowed edge internally. Thus these are certificates of upper bounds
on cardinality, without assuming any coloring heuristic is optimal.

It suffices to check all distinct values `rho = U_j <= 10`, including 10.
Between consecutive thresholds, the allowed graph is constant up to its
upper endpoint, and density increases with rho. With `rho = t/m`, each
threshold yields the exact bound

```
delta² <= (c+1)² * (t/(4m))^55 / (29^7 m).
```

Every bound in all 24 families is strictly smaller than
`(5/2)^55 / (2*3^24)`. The largest is the 103,140.241156608598… ceiling
in case 12. All comparisons use rational squares; displayed square roots
are decimal approximations. [results.json](results.json) contains the
full per-family bounds and thresholds attaining them.

## Reproduction

From the repository root, with NumPy and python-flint installed:

```sh
python research/mw56_projection55_20261116/verify.py
```

This repeats the parent reconstruction, all section determinants, all
distance witnesses, every coloring, all density comparisons, and the five
explicit rejections. It checks the stored result rather than rewriting it.

To rerun the initial discovery screen, additionally install fpylll and run:

```sh
python research/mw56_projection55_20261116/search.py --output /tmp/mw56-screen.json --limit 24
```

Choose a fresh output path. This script reproduces the initial reduction
and short-vector walk; it does not rerun the later targeted enumerations
that supplied additional certificate vectors. Its output is a search
report, not a packing certificate. The archived final witnesses remain
independently verifiable without repeating their discovery.
