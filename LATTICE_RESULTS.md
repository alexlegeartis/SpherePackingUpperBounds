# Lattice-packing density upper bounds

These bounds apply **only to lattice packings**. They do not replace the unrestricted sphere-packing results in [RESULTS.md](RESULTS.md). General-packing upper bounds apply to lattices too; this table does not relabel those results as new lattice-specific improvements.

| Dimension | Known lattice lower L | Cohn snapshot upper U | New certified upper U′ | Upper decrease | Gap decrease | Method |
|---:|---:|---:|---:|---:|---:|---|
| 14 | ≈0.0216240960824471 | 0.0352392079482642 | **≤0.034712755778** | 1.4939% | 3.8667% | Near-frame arithmetic + mod-2 class budgets ([proof/replay](verifications/lattice/dim14/README.md)) |
| 15 | ≈0.0168575706567626 | 0.0256118266556343 | **≤0.025511135431** | 0.3931% | 1.1502% | Near-frame arithmetic + mod-2 class budgets ([proof/replay](verifications/lattice/dim15/README.md)) |
| 16 | ≈0.0147081643974308 | 0.0187461205770637 | **≤0.018655117428** | 0.4855% | 2.2537% | Near-frame arithmetic + mod-2 class budgets ([proof/replay](verifications/lattice/dim16/README.md)) |
| 17 | ≈0.0088113191823211 | 0.0138839312698563 | **≤0.013832417384** | 0.3710% | 1.0155% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim17/README.md)) |
| 18 | ≈0.0059283687184694 | 0.0103704833569080 | **≤0.010343803448** | 0.2573% | 0.6006% | Near-frame arithmetic + mod-2 class budgets ([proof/replay](verifications/lattice/dim18/README.md)) |
| 19 | ≈0.0041208062797686 | 0.0078244736036902 | **≤0.007808430327** | 0.2050% | 0.4332% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim19/README.md)) |
| 25 | ≈0.0006772120097731 | 0.0013293080953711 | **≤0.001325232433** | 0.3066% | 0.6250% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim25/README.md)) |
| 26 | ≈0.0002692200504338 | 0.0009151191756312 | **≤0.000907209997** | 0.8643% | 1.2245% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim26/README.md)) |
| 27 | ≈0.0001286752812026 | 0.0006313421536495 | **≤0.000622938953** | 1.3310% | 1.6717% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim27/README.md)) |
| 28 | ≈0.0000697587366165 | 0.0004370005973032 | **≤0.000429639101** | 1.6846% | 2.0045% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim28/README.md)) |
| 29 | ≈0.0000278789874568 | 0.0003038425430867 | **≤0.000298025039** | 1.9146% | 2.1081% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim29/README.md)) |
| 30 | ≈0.0000144286438425 | 0.0002124025969562 | **≤0.000207932946** | 2.1043% | 2.2577% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim30/README.md)) |
| 31 | ≈0.0000118377651859 | 0.0001494085876507 | **≤0.000146087217** | 2.2230% | 2.4143% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim31/README.md)) |
| 32 | ≈0.0000110407493088 | 0.0001057644820174 | **≤0.000104134823** | 1.5408% | 1.7204% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim32/README.md)) |
| 33 | ≈0.0000041406882896 | 0.0000755167889394 | **≤0.000073827123** | 2.2375% | 2.3673% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim33/README.md)) |
| 34 | ≈0.0000017669738891 | 0.0000543081490854 | **≤0.000053148401** | 2.1355% | 2.2073% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim34/README.md)) |
| 35 | ≈9.461904151E-7 | 0.0000393266736309 | **≤0.000038509364** | 2.0783% | 2.1295% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim35/README.md)) |
| 36 | ≈7.330347783E-7 | 0.0000286680608768 | **≤0.000028235617** | 1.5085% | 1.5480% | Global cross-class coupling + mod-2 class budgets ([proof/replay](verifications/lattice/dim36/README.md)) |

## Comparison scope and sources

The comparison columns reproduce the **lattice** section of a user-supplied snapshot of [Henry Cohn’s sphere-packing table](https://cohn.mit.edu/sphere-packing/), received on 2026-11-16. The snapshot contains no extraction date; this is not a claim that the live page was fetched or is unchanged. Every new upper bound here is smaller than the corresponding upper in that snapshot. Its listed previous upper bounds in these dimensions are attributed to [Nina Zubrilina, *Improved lattice-packing bounds in small dimensions*](https://arxiv.org/abs/2610.35830).

The snapshot gives more digits than the outward claims in the original source certificates. The JSON retains both the snapshot comparison upper and the original certificate-claim upper. Original comparison certificates are preserved from the [pinned upstream repository](https://github.com/AstroNinja389/lattice-packing-bound/tree/340a0e93bb95942c2b1b80092be43149504970eb/certificates). The 14D, 15D and 18D prior certificates are exactly replayed because those arithmetic proofs require them; later global-coupling results have no prior-density assumption.

Known lattice lower values are the displayed values in the same snapshot, not the generally larger nonlattice entries in its unrestricted table. They are rounded/truncated comparisons, not exact decimal identities or independently exhaustive current-record claims. In 36D the snapshot already cites this repository’s lattice construction. Gap decrease is `(U−U′)/(U−L)` using these displayed values; upper decrease is `(U−U′)/U`. These are comparison statistics, not new theorems.

This establishes improvements over the supplied snapshot and the named source benchmarks, not optimality or an exhaustive worldwide-priority guarantee. The snapshot credits this repository’s general-packing upper bounds in 23D and 37D–48D; those remain in the separate unrestricted table rather than being duplicated here.

### Results deliberately not included as current records

- 11D–13D calculations were superseded by the newer claimed computer-assisted results listed in the snapshot. See [literature status](verifications/lattice/LITERATURE_STATUS.md). Their complete separate proof chains have not been independently replayed here.
- Earlier tiny refinements in 20D, 26D and 32D were not published because those old proof packages were not recovered. The new 26D and 32D certificates above are different, independently audited global-coupling results; 20D remains excluded here pending a complete published package.
- 9D and 10D exploratory values were weaker than known lattice bounds.

## Method and reproducibility

The 14D, 15D, 16D and 18D certificates combine Zubrilina’s modulo-two class framework with arithmetic near-frame constraints, exact torsion exclusions, and exhaustive first-n profile dual bounds. The determinant step overlaps [Géraud-Stewart and Naccache, Appendix D, Lemma 32](https://eprint.iacr.asia/2026/2826.pdf); it is attributed and not claimed as new.

The 17D, 19D and 25D–36D certificates use global short-class angular coupling, continuous Gegenbauer polynomial inequalities, and exact projective multiplicity bounds. They have no arithmetic trace-gap or prior-density assumption. Proof documents explain the continuous geometry and source-count normalization. Every published result has a complete executable exact proof.

Run the separate lattice suite without changing the historical unrestricted suite:

```sh
python verifications/lattice/replay.py
# Or one dimension:
python verifications/lattice/replay.py --dimension 36
```

Python 3.10+ is required; the earlier dimensions use SymPy. The default 25D–36D full verifiers use only the standard library, with exact compactified-half-line Bernstein sign proofs. Separate optional SymPy verifiers are included. Public `python-flint` is an optional speed improvement, not a private dependency.

Large JSON inputs use deterministic gzip. The replay checks stored and decompressed SHA256 hashes, materializes a temporary copy, and runs the complete proof there. Audited certificate bytes are preserved. Do not use `python -O`. Hosted Actions have not been enabled or dispatched.
