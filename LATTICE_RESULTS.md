# Lattice-packing density upper bounds

These bounds apply **only to lattice packings**. They do not replace the unrestricted sphere-packing results in [RESULTS.md](RESULTS.md). General-packing upper bounds apply to lattices too; see that separate table rather than treating its lattice applicability as a new result.

| Dimension | Named lattice lower comparison L | Earlier lattice upper U | New certified upper U′ | Upper decrease | Gap decrease | Certificate |
|---:|---:|---:|---:|---:|---:|---|
| 14 | ≈0.0216240960824471 | 0.0352392080 | **≤0.034712755778** | 1.4939% | 3.8667% | [Proof and replay](verifications/lattice/dim14/README.md) |
| 15 | ≈0.0168575706567627 | 0.0256118267 | **≤0.025511135431** | 0.3931% | 1.1502% | [Proof and replay](verifications/lattice/dim15/README.md) |
| 16 | ≈0.0147081643974308 | 0.0187461205771 | **≤0.018655117428** | 0.4855% | 2.2537% | [Proof and replay](verifications/lattice/dim16/README.md) |
| 18 | ≈0.00592836871846942 | 0.0103704834 | **≤0.010343803448** | 0.2573% | 0.6006% | [Proof and replay](verifications/lattice/dim18/README.md) |

## Comparison scope and sources

Earlier upper bounds are from [Nina Zubrilina, *Improved lattice-packing bounds in small dimensions*](https://arxiv.org/abs/2610.35830) and its [certificate repository](https://github.com/AstroNinja389/lattice-packing-bound/tree/340a0e93bb95942c2b1b80092be43149504970eb/certificates). The displayed 14D, 15D and 18D values are outward claims exactly replayed from those source certificates; 16D uses the earlier certificate value cited in the research comparison. The lower comparisons are explicit laminated/Barnes–Wall lattices from the [Nebe–Sloane lattice catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/density.html): center densities 1/(16√3), 1/(16√2), 1/16 and 1/(8√3), respectively. Displayed packing fractions are rounded; these are verified known lattice lower bounds, not an exhaustive current lower-bound priority claim. Percentage gap reduction is `(U−U′)/(U−L)` using the displayed L; it is a comparison statistic, not a new theorem.

These are verified improvements over the named benchmarks. They are not optimality results or an exhaustive worldwide-priority guarantee. The 16D result also improves the earlier reported refinement 0.018746011555.

### Results deliberately not included as current records

- 11D–13D calculations were superseded as record comparisons by newer claimed computer-assisted results. See [literature status](verifications/lattice/LITERATURE_STATUS.md). The authors report full proofs; this repository has not independently replayed their complete separate proof chain.
- Earlier reported refinements in 20D, 26D and 32D are not included here because their full proof packages have not been recovered and replayed for this publication.
- 9D and 10D exploratory values were weaker than known lattice bounds and are not records.

## Method and reproducibility

The new certificates combine Zubrilina’s modulo-two class framework with arithmetic near-frame constraints, exact torsion exclusions, and exhaustive first-n radius-profile dual bounds. The determinant step overlaps [Géraud-Stewart and Naccache, Appendix D, Lemma 32](https://eprint.iacr.asia/2026/2826.pdf); it is attributed in the proof documents and is not claimed as a new determinant lemma.

Run the separate lattice suite, without modifying the historical unrestricted suite:

```sh
python verifications/lattice/replay.py
# Or one dimension:
python verifications/lattice/replay.py --dimension 16
```

Requires Python 3.10+ and SymPy. Large JSON inputs are stored as deterministic gzip files. The replay materializes a temporary copy, verifies every decompressed SHA256, and runs the exact proof programs there. Audited certificate bytes are unchanged. Do not use `python -O`. Hosted Actions have not been enabled or dispatched.
