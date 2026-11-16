# 41D lattice-only bound

Every 41-dimensional lattice packing has density **<0.000004531256243**. This is not a claim for unrestricted packings.

Reproduce with Python 3.12+, the pinned `requirements.txt`, and:

```sh
python verifications/lattice/replay.py --dimension 41
```

All compressed and decompressed file hashes are checked before replay in a temporary directory. Assertions are required. No optimizer, network access, density-bound premise or classification assumption is needed.

Certificate SHA256: `087f3bae78580cfc01aa26c6e3ff2d0e31038bf095b6806121172b9ef35c1773`.

The proof combines a cyclic PSD kernel and forced collinear triples with mod-2 class profiles, conditioned Gram-moment capacities and angular pair-sum cuts. See [class/AA foundations](CLASS_PROFILES.md), [kernel bridge](CYCLIC_CLASS_BRIDGE.md), [collinear proof](COLLINEAR_CYCLIC_CREDIT.md), [all-short constraints](ALL_SHORT_PROOF.md), [angular inequalities](ANGULAR_PROOF.md) and [projective caps](PROJECTIVE_CAP_PROOF.md). Historical 38D illustrations in component notes are not assumptions for this dimension; all saved coefficients and finite-grid constraints are replayed.

Two independent radial/profile implementations recheck Fourier reconstruction, both global signs, every one of 240 mixed bins, all 4 AAcuts, 193 class branches, weighted inequalities and normalization. Their envelope algorithms differ, but both use Arb. Saved audit reports are not trusted proof inputs.

This uses a NEW dimension 41 angular polynomial (eta=.30 with margin .01), SHA 53033e2d2b55008c352a252bf2ee4f3f6c41ad10a71a87fab65ae9a3b2a49519, not the older unrestricted polynomial. Two complete angular implementations are replayed: the original primary rectangle gate and the independent centered-gradient gate. See the [centered-gradient derivation](MEANVALUE_CYCLIC_GATE.md). The primary replay takes several minutes. Both gates include exact harmonic coefficients, cap/transition, all-scale and infinite-tail checks. Their interval backend is shared Arb.


The profile is h(t)=t^96 exp(-120(t-1)), and the origin-block ratio is exactly **7/10**. The [origin-ratio proof](ORIGIN_RATIO_PROOF.md) supersedes the ratio-one formulas in the generic bridge notes: the spatial negative cross term is multiplied by r, while the objective origin cost is r²C. All these terms and the combined tail are checked exactly. The underlying mod-2 framework extends [Zubrilina (2026)](https://arxiv.org/abs/2610.35830); source cyclic methodology and generic angular proofs build on Alex Kravatskiy's prior repository packages. The preceding unrestricted density bound is only a comparison, not a theorem premise. Lattice lower/gap cells remain blank pending separately approved source updates.

Endpoint convention: radial bins are right-open [E_j,E_{j+1}) and cumulative cutoffs are strict. Closed intervals are used only for conservative envelope checks. Thus A=0 at short_upper+long_upper≤4 is valid. Vectors at the final squared norm9/4 belong to the nonpositive tail. The j=−1 branch is the no-short-vector class branch, not a literally empty class.
