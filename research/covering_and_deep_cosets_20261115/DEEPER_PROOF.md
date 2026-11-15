# Improved explicit lattices in 49–52 dimensions

The old four constructions are superseded by the explicit lattices here.
The new ingredient is one singular P48q half-coset of squared depth at
least **4**, combined with shallower cosets and a different vertical
metric. All four lattices have minimum exactly **6**.

## Horizontal certificate

The [integer certificate](explicit_p48/certificate.npz) gives four binary
residue rows B in the original P48q basis and in the same reduced basis
as the preceding package. The exact setup audit binds the original
Gram to [the catalogue snapshot](../range48_96_20261115/explicit_p48/P48q.html),
checks the unimodular basis change, and verifies that all 15 nonzero
binary combinations represent distinct nonzero classes.

For coefficients `(a0,a1,a2,a3)` the parent quadratic type is

```
q(a) = a1+a2+a1*a2 mod 2.
```

The special class `a=(1,0,0,0)` has no representative of norm at most
12: a **complete** congruence-class enumeration verifies this. Its
representative norms are divisible by 4, so its half-coset has squared
minimum at least `16/4=4`.

For each of the other fourteen classes, a complete enumeration rules
out norm at most 10. Hence their half-coset minima are at least 3 for
singular classes and 7/2 for nonsingular classes. The published P48q
minimum-6 theorem is inherited, not re-enumerated.

All fifteen enumerations completed in the release replay, visiting
**638,998,811 nodes**. A deliberately short class was correctly rejected
with an exact short-vector witness. The full input files, reports and
compiler/source hashes are in [the replay archive](explicit_p48/replay/manifest.json).

## Vertical metrics

Use the leading r-by-r block M_r, r=1,2,3,4, of

```
M4 = [[2,  1,   1,    0   ],
      [1,  5/2, 5/4, -3/2 ],
      [1,  5/4, 5/2, -3/2 ],
      [0, -3/2,-3/2,  3   ]].
```

For every nonzero integer coefficient vector z, this metric satisfies

| Binary class | Required vertical squared norm |
| --- | ---: |
| zero | at least 6 |
| special class `(1,0,0,0)` | at least 2 |
| other singular class | at least 3 |
| nonsingular class | at least 5/2 |

The checker proves these inequalities by complete small-dimensional
integer enumeration. This is exhaustive: every vector with norm below
6 obeys `z_i^2 < 6*(M_r^-1)_ii`. An exact rational bound therefore
confines all possible violations to the checked coefficient box.
Positive definiteness and every inverse and determinant are exact.

The metric was found by fitting the differing half-coset depths. In
three dimensions, the candidate has determinant 55/8. A weighted
Voronoi search then supplied the rational fourth-layer offset
`w=(6,-6,-6)/11` and height square `15/11`. These give M4 and
determinant 75/8. The proof uses the exact final metric and the complete
integer audit, not an assertion that the numerical search was globally
optimal.

## Full lattices and density

Let B_r consist of the first r original residue rows. Generate the
lattice by the parent L and `(b_i/2,v_i)`, where v_i have Gram M_r.
Every nonzero layer has horizontal plus vertical norm at least 6:
the three nonzero binary cases sum to `4+2`, `3+3`, or `7/2+5/2`.
All-even nonzero layers meet the vertical bound on their own. The
parent attains 6, so the complete lattice minimum is exactly 6.

The full Gram is

```
[[ G,           G B_r^T/2              ],
 [ B_r G/2,     B_r G B_r^T/4 + M_r    ]].
```

The auditor reconstructs this matrix and checks its determinant
against det(M_r). Thus `delta_(48+r)^2=(3/2)^(48+r)/det(M_r)`.

| Dimension | Gram determinant | Center density |
| --- | ---: | ---: |
| 49 | 2 | 14,578.7688119179… |
| 50 | 4 | 12,625.5841470212… |
| 51 | 55/8 | 11,794.8047134541… |
| 52 | 75/8 | 12,370.4955459097… |

The [JSON result](deeper_results.json) includes the actual residue rows,
vertical Grams and exact rational density squares. Each is strictly
stronger than the corresponding checkpoint construction.

## Reproduce

Quick exact setup audit and checks of the saved completed reports:

```sh
python research/covering_and_deep_cosets_20261115/verify_deeper.py
```

Fresh full replay of all 15 cosets and the negative control:

```sh
python research/covering_and_deep_cosets_20261115/verify_deeper.py --output /tmp/deeper_p48_replay
```

The directory must be new. This needs NumPy, SymPy, python-flint and
GNU C++. The enumerator adapts the preceding audited P48 implementation,
compiled separately with cutoffs 10 and 12. It uses outward-rounded
IEEE-754 binary64 inequalities, exact dyadic LDL enclosures and exact
integer leaf checks, with `-fno-fast-math -ffp-contract=off` and range
guards. These are the same documented arithmetic/compiler assumptions
as the preceding P48 certificate. The quick audit does not rerun the
enumerations; the `--output` command does.
