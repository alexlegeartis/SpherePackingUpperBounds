# Stronger lower bounds in thirteen dimensions

**Six larger gains and seven smaller covering refinements**, all within 48D–96D. The new explicit 53D lattice replaces a conservative section estimate. The other results strengthen our existing dated best-known bounds. The [current source review](../covering_and_deep_cosets_20261115/PRIORITY.md) now registers 39 construction dimensions; no stronger published lower bound was located in that review.

All densities below are approximate **center densities**. The incremental gains compare with the exact bounds pinned at commit `22bda7e` in [checkpoint.json](checkpoint.json). The main repository table separately retains the named external comparisons.

| Dimension | Previous bound | New bound | Increase over preceding bound | Proof |
| --- | ---: | ---: | ---: | --- |
| 50 | 12625.58414702117 | **14578.768811917893** | **15.470054%** | [Deep plane + fitted metric; explicit minimum-6 lattice](P48_PROOF.md) |
| 51 | 11794.80471345414 | **13186.992560309189** | **11.803399%** | [Deep plane + fitted metric; explicit minimum-6 lattice](P48_PROOF.md) |
| 52 | 12370.49554590971 | **13773.364524023099** | **11.340443%** | [Deep plane + fitted metric; explicit minimum-6 lattice](P48_PROOF.md) |
| 53 | 17624.6617556196 | **20890.804898853429** | **18.531664%** | [Primitive dual normal; explicit minimum-10 lattice](SECTION53_PROOF.md) |
| 77 | 61083979320.88889 | **70533703808.177467** | **15.470054%** | [Adaptive exceptional coset above D5; lattice existence proof](ADAPTIVE_PROOF.md) |
| 85 | 1453530529565.988 | **1598670440672.3403** | **9.985336%** | [Adaptive exceptional coset + fitted metric; lattice existence proof](ADAPTIVE_PROOF.md) |

The 50D–53D constructions are explicit lattices. The 77D/85D constructions prove existence of lattices; horizontal representatives remain unknown.

The [subsequent MW92 continuation](../mw92_sections_20261115/) adds six new construction dimensions, 86D–91D.

## Comparisons with named earlier work

| Dimension | External comparison | Increase over that comparison |
| --- | --- | ---: |
| 50 | P48q × A2, **derived** from Bos–Conway–Sloane (1982) §5 | **200.000000%** |
| 51 | P48q × A3, **derived** from Bos–Conway–Sloane (1982) §5 | **343.129368%** |
| 52 | Two Blichfeldt sections of Leterrier54, **derived** from Leterrier (2022) Corollary 1.4 | **162.973303%** |
| 53 | One Blichfeldt section of Leterrier54, **derived** from Leterrier (2022) Corollary 1.4 | **19.152528%** |
| 77 | Bare L80 section with audited rank 3 normal, **derived** from Bachoc–Nebe (1998) | **214.269681%** |
| 85 | BCS88 perpendicular to √2A3, **derived** from Bos–Conway–Sloane (1982) Theorems 1–2 | **45.398230%** |

[Source links and exact comparison formulas](../covering_and_deep_cosets_20261115/PRIORITY.md) · [Machine-readable external comparisons](../covering_and_deep_cosets_20261115/literature_comparisons.json). Derived numbers are not attributed as printed records to the source authors.

## The two main construction improvements

**Preserve a whole deep plane.** The P48q search found three compatible depth-4 half-cosets, then twelve compatible nonsingular classes of depth at least 7/2. Fitting the vertical Gram to these depths gives the three explicit improvements in 50D–52D. All fifteen coset exclusions completed, totaling **1,002,136,588 nodes**. A fresh full rerun of the same audited enumerator completed with matching node counts and passed the negative control. Both report sets are retained; the rerun is not an independent enumeration implementation.

**Adapt to one exception instead of requiring none.** A weighted integer bad-coset cost has expectation below 2. Hence one can select a model with at most one shallow class, and choose the vertical lattice to absorb that class. In 77D, move it to a longer D5 parity class. In 85D, use a fitted rational metric and D5 root isometries. This avoids throwing away an otherwise useful subspace when the old zero-exception union bound fails. [Complete argument](ADAPTIVE_PROOF.md).

The separate 53D construction uses an explicit primitive dual normal of squared norm 10/3 in Leterrier's displayed 54D Gram. HNF constructs the complete kernel, its determinant is 10·3^24, and a norm 10 vector proves minimum attainment. [Source and exact proof](SECTION53_PROOF.md).

## Smaller covering refinements

A global cubic Taylor lower polynomial improves the finite Coxeter–Few–Rogers simplex estimate. These seven gains are modest and are listed separately. [Proof](COVERING_PROOF.md).

| Dimension | New center density | Increase over preceding bound |
| --- | ---: | ---: |
| 57 | 163093.23475965115 | 0.022355% |
| 81 | 1066394783124.8914 | 0.013686% |
| 82 | 1040341684972.145 | 0.025997% |
| 93 | 8743033982260848.9 | 0.011284% |
| 94 | 8584803057662693.3 | 0.020953% |
| 95 | 8868241996705554.8 | 0.030350% |
| 96 | 10395177038683104.0 | 0.039348% |

The updated 96D bound is **118.640977% above the explicitly reported BCS1982 construction** `3^48/2^24`.

## Verification and reproducibility

```sh
python research/adaptive_layers_20261115/verify_all.py
python research/adaptive_layers_20261115/verify_p48.py --replay /tmp/new_p48_replay
python research/adaptive_layers_20261115/verify_section53.py --source-pdf /path/to/EPFL_TH10417.pdf
```

The first command audits exact arithmetic, complete saved enumeration reports and the second full P48 report set. The second recompiles and reruns the complete enumeration into a new directory. NumPy and python-flint are required; the inherited P48 setup also uses SymPy. Do not run with `-O`. The optional metric discovery script uses SciPy, which is not needed by any proof auditor.

Published parent minima, the computational 54D Gram, theta-series and covering theorems are stated external inputs. The P48 enumeration uses the documented outward-rounded IEEE-754/compiler assumptions. All other new matrix and finite layer checks use exact integer/rational arithmetic. [Sources and hashes](sources.json) · [Search budgets](discovery_summary.json) · [All results](results.json).
