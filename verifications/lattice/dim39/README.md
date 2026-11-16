# 39D lattice-only bound

Every 39-dimensional lattice packing has density **<0.000009086501484**. This is not an unrestricted-packing claim.

The proof combines a cyclic PSD kernel, forced collinear triples, mod-2 class profiles, conditioned Gram-moment inequalities and angular pair-sum cuts. It has no prior density bound or quotient-classification premise. No worldwide priority claim is made.

## Replay

Install Python 3.12+ and the pinned packages in `requirements.txt`, then run from the repository root:

```sh
python verifications/lattice/replay.py --dimension 39
```

The manifest checks compressed/decompressed SHA256 values and materializes a temporary copy. `run_all.py` replays both original cyclic gates, the first independent radial/profile audit, and the second independent complete audit. The latter rechecks all weighted-class constraints at every used cutoff. Assertions are mandatory; `python -O` is rejected. No optimizer or network is needed.

Certificate SHA256: `32915e94df222f092af74294917a4896df8a3f85d6c20db17e1ba085e368bc31`. Both independent implementations check all240 bins, 2 angular cuts and193 class branches, Fourier reconstruction, global polynomial signs, complete mixed-exponential envelopes, profile dual domination and final normalization. They use distinct envelope algorithms but share the Arb backend. Saved reports are provenance, not trusted proof inputs.

Normalization is an exact squared comparison using a rational Machin lower bound for pi (20 and21 alternating-series terms).

## Proof and attribution

See [foundations](CLASS_PROFILES.md), [collinear residual](COLLINEAR_CYCLIC_CREDIT.md), [kernel bridge](CYCLIC_CLASS_BRIDGE.md), [all-short constraints](ALL_SHORT_PROOF.md), [angular coupling](ANGULAR_PROOF.md) and [projective caps](PROJECTIVE_CAP_PROOF.md). The formulas are dimension-generic. Historical38D numerical illustrations in component notes are not additional hypotheses; this dimension's coefficients and all finite-grid checks are replayed anew.

The identical dimension-specific angular polynomial and profile h(t)=t^96 exp(-120(t-1)) come from Alex Kravatskiy's [original unrestricted package](../../dim39/refinement_profile20261114/), with original shared gate scripts pinned to commit107d86d8d1c5c4e029ef10845b22519738592620. `cyclic/SOURCE_MANIFEST.json` records original hashes. Its metadata-only certificate is not an old radial certificate. The generic angular proof text was copied from the same commit's dim37 lattice package. The mod-2 class framework extends [Nina Zubrilina (2026)](https://arxiv.org/abs/2610.35830).

The comparison predecessor is an unrestricted upper bound, hence also applies to lattices. It is used only for percentage comparisons, not assumed by this proof. Lattice lower and gap cells are omitted pending a verified lattice-only lower source.
