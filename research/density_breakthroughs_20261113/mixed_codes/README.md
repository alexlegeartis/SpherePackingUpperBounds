# Mixed quaternary/ternary hexagonal constructions

**No new packing was obtained.** One proposed34D route is excluded exactly;
three rank-three parity searches failed to find a feasible construction.

The local squared-distance costs are those in the existing 36D provenance:
0 for `(0,0)`,4 for `(0,nonzero ternary)`,3 for `(nonzero quaternary,0)`,
and 1 when both are nonzero. For a GF(4) outer code of length n and size M,
and a rank-r ternary parity matrix with no zero column, the proposed
minimum-eight construction has center density

```
M / 3^(n/2+r).
```

When the outer minimum is 6 or 7, every support S of a weight6/7 word must avoid
a full-support signed dependence among the ternary columns:
`sum_{i in S} epsilon_i H_i != 0` for every `epsilon_i in {+1,-1}`.
The searches enforce this necessary-and-sufficient obstruction to distances
below 8 for these outer codes. A zero parity column is inadmissible because
it allows a ternary weight-one word with local squared cost4.

## Exact obstruction for the length17, dimension 9 outer code

The archived Grassl construction punctures the known `[18,9,8]` GF(4) code.
Two independent complete enumerations found 262,144 words, minimum 7,
1224 words of weight7, and 408 distinct weight7 supports. Every coordinate
pair belongs to exactly 63 of these supports: a `2-(17,7,63)` design.

A nonzero column in F3^2 has one of four projective colors. Exact enumeration
of all seven-column multiplicity patterns shows that a signed dependence
is avoided only when the four color counts are a permutation of `(6,1,0,0)`.
Each support must therefore contain15 monochromatic coordinate pairs.
Counting these pairs over the design would give

```
63 * (total monochromatic coordinate pairs) = 408*15,
total monochromatic coordinate pairs = 680/7,
```

which is impossible. This excludes **all nonzero rank-at-most-two ternary
parity assignments for this particular outer code**, not arbitrary34D
packings or every quaternary code.

`verify.py` reconstructs all 262,144 words, the support design, and the
projective-color classification. `../sections/audit_n17.py` is an independent
implementation using polynomial multiplication in GF(4).

## Rank-three searches

The `[20,13,6]` generator is from the archived Grassl/Edel–Bierbrauer source.
Shortening yields the length19/dimension 12 and length18/dimension 11 inputs.
`enumerate.cpp` visits every binary Gray-code word and checks the exact weight
enumerators, including uniqueness of the zero word.

| Outer code | Exact word count | Distinct weight6/7 supports | Best remaining violated supports |
|---|---:|---:|---:|
| [18,11,6] | 4,194,304 | 1,757 | 384 |
| [19,12,6] | 16,777,216 | 2,716 | 619 |
| [20,13,6] | 67,108,864 | 4,100 | 979 |

These remaining counts are independently recomputed using exact signed sums.
The rank-three searches are heuristic; their failures prove no nonexistence
result. In particular, the attractive formal densities must not be reported
as packing lower bounds.

## Reproduction

From the repository root, run `build_tools.sh` in the parent directory, then:

```bash
R=research/density_breakthroughs_20261113
python "$R/mixed_codes/prepare.py"
for name in n17k9 n18k11 n19k12 n20k13; do
  "$R/bin/mixed_enumerate" "$R/mixed_codes/$name.input" "$R/mixed_codes/$name.supports"
done
"$R/bin/mixed_search" "$R/mixed_codes/n18k11.supports" 3 90 4103 "$R/mixed_codes/n18k11_parity.json"
"$R/bin/mixed_search" "$R/mixed_codes/n19k12.supports" 3 45 4019 "$R/mixed_codes/n19k12_parity.json"
"$R/bin/mixed_search" "$R/mixed_codes/n20k13.supports" 3 45 4020 "$R/mixed_codes/n20k13_parity.json"
python "$R/mixed_codes/verify.py"
```

The finite enumerations are deterministic. Time-budgeted search outputs may
vary with machine speed; saved outputs can be checked without rerunning them.
