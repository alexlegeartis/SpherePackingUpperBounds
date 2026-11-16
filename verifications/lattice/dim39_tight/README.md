# 39D lattice-only bound

Every 39-dimensional lattice packing has density **<0.000009068287259**. This is not a claim for unrestricted packings.

Reproduce with Python 3.12+, the pinned `requirements.txt`, and:

```sh
python verifications/lattice/replay.py --dimension 39
```

All compressed and decompressed file hashes are checked before replay in a temporary directory. Assertions are required. No optimizer, network access, density-bound premise or classification assumption is needed.

Certificate SHA256: `34966e32a29f93207c36e1b86018ab0b367a9c56cc5b5695059a6cba3bea8b9e`.

The proof combines a cyclic PSD kernel and forced collinear triples with mod-2 class profiles, conditioned Gram-moment capacities and angular pair-sum cuts. See [class/AA foundations](CLASS_PROFILES.md), [kernel bridge](CYCLIC_CLASS_BRIDGE.md), [collinear proof](COLLINEAR_CYCLIC_CREDIT.md), [all-short constraints](ALL_SHORT_PROOF.md), [angular inequalities](ANGULAR_PROOF.md) and [projective caps](PROJECTIVE_CAP_PROOF.md). Historical 38D illustrations in component notes are not assumptions for this dimension; all saved coefficients and finite-grid constraints are replayed.

Two independent radial/profile implementations recheck Fourier reconstruction, both global signs, every one of 240 mixed bins, all 2 AAcuts, 193 class branches, weighted inequalities and normalization. Their envelope algorithms differ, but both use Arb. Saved audit reports are not trusted proof inputs.

This is a tight-envelope refinement of the preserved [previous package](../dim39/README.md): kernel, radial polynomials and theorem assumptions are unchanged, but every gate is replayed. The unchanged dimension-specific kernel and gates originate from the pinned repository commit 107d86d8d1c5c4e029ef10845b22519738592620.

The profile is h(t)=t^96 exp(-120(t-1)), and the standard origin-block ratio is exactly 1. No generalized-origin-block result is claimed. The underlying mod-2 framework extends [Zubrilina (2026)](https://arxiv.org/abs/2610.35830); source cyclic methodology and generic angular proofs build on Alex Kravatskiy's prior repository packages. The preceding unrestricted density bound is only a comparison, not a theorem premise. Lattice lower/gap cells remain blank pending separately approved source updates.

Endpoint convention: radial bins are right-open [E_j,E_{j+1}) and cumulative cutoffs are strict. Closed intervals are used only for conservative envelope checks. Thus A=0 at short_upper+long_upper≤4 is valid. Vectors at the final squared norm9/4 belong to the nonpositive tail. The j=−1 branch is the no-short-vector class branch, not a literally empty class.
