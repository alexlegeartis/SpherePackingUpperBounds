# 38D lattice-only bound

Every 38-dimensional lattice packing has density **< 0.000013058679142**.
This is not a bound for arbitrary packings. The method combines a cyclic PSD kernel with forced collinear triples, mod-2 class profiles, conditioned Gram-moment capacities, and pair-sum angular inequalities.

## Replay

Install Python 3.12+, `python-flint` and `sympy`. From the repository root run:

```sh
python verifications/lattice/replay.py --dimension 38
```

The materializer checks compressed and decompressed SHA256 values and runs in a temporary directory. Two pinned cyclic angular gates, class inequalities, and both independent complete certificate implementations must pass. No optimizer, network access or external density/classification premise is needed. Assertions are mandatory; optimized Python mode is rejected.

Frozen certificate SHA256: `34888f6e673b2369bf69cd506e6fab0c3cf8612a03c3617ff85ca0d1678a9eda`.

The two independent implementations check exact Fourier reconstruction, both global polynomial rays, all 240 mixed radial bins, four angular cuts, 193 class branches, and rational even-dimensional density normalization. Saved audit reports are provenance, not substitutes for replay.

## Sources and proof dependencies

The angular polynomial and all-scale profile are copied unchanged from Alex Kravatskiy's [unrestricted 38D package](../../dim38/refinement_profile20261114/) and [shared gates](../../profile_enrichment_20261114/) at commit `107d86d8d1c5c4e029ef10845b22519738592620`. Their original hashes appear in `cyclic/REPLAY_MANIFEST.json`. The small cyclic certificate contains only dimension/profile metadata, not an old radial proof. The mod-2 class framework extends [Nina Zubrilina (2026)](https://arxiv.org/abs/2610.35830).

See [method inventory](METHOD.md), [kernel bridge](CYCLIC_CLASS_BRIDGE.md), [collinear residual](COLLINEAR_CYCLIC_CREDIT.md), [independent derivation](COLLINEAR_RESIDUAL.md), and [all-short inequality](ALL_SHORT_PROOF.md), [class capacities and weighted ranks](CLASS_PROFILES.md), [angular coupling](ANGULAR_PROOF.md), and [projective cap proof](PROJECTIVE_CAP_PROOF.md). The last two generic proofs are copied unchanged from the pinned repository's dim37 package. Historical statements in these component notes that a density certificate remained necessary are discharged by this package's two full replays.

The comparison upper 0.0000145182018536 is the preceding unrestricted bound, also valid for lattices, and is used only for the 10.0531% comparison, not as a theorem premise. No worldwide priority claim is made. The lower/gap comparison is omitted pending a verified lattice-only lower source.

## Tight-envelope refinement

This package tightens only the radial envelope certificates and profile duals of the [preceding 38D package](../dim38/README.md), which remains immutable. The kernel, radial polynomials, angular cuts and theorem assumptions are unchanged. Both independent implementations and all source gates are nevertheless replayed in full. The exact input equalities are documented in `INHERITANCE.json`; they are provenance, not a substitute for replay.

Endpoint convention: radial bins are right-open [E_j,E_{j+1}) and cumulative cutoffs are strict. Closed intervals are used only for conservative envelope checks. Thus A=0 at short_upper+long_upper≤4 is valid. Vectors at the final squared norm9/4 belong to the nonpositive tail. The j=−1 branch is the no-short-vector class branch, not a literally empty class.
