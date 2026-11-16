# 44D lattice-only bound

Every 44-dimensional lattice packing has density **<0.000001783944991**. This is not an unrestricted-packing claim.

Install Python 3.12+ and the pinned requirements. From the repository root run `python verifications/lattice/replay.py --dimension 44`; alternatively run `python run_all.py` inside a materialized package. No optimizer, network or prior density/classification assumption is needed. Assertions must remain enabled.

Certificate SHA256: `2b373f4af6df2b750836328262f79af78f2f72ef7ccfd00e659913c3b3ad923d`.

The complete centered-gradient angular gate is shared by the two independent radial/profile audits. Its analytic derivative derivation was independently reviewed; this is not a claim of two independent angular implementations. The runner replays the whole gate, then both independent radial/profile implementations, including all 240 bins, 1 AA cut, 193 class branches, Fourier/global signs, weighted rank inequalities, exact envelopes and normalization. Both algorithms use Arb.

The new 44D angular polynomial has SHA b79ca676bdbcc8cb5739961d4c46bbadeda662e8042a8b2a0d27df6777c2483a and its complete cyclic domain is rechecked.

See [class foundations](CLASS_PROFILES.md), [angular proof](ANGULAR_PROOF.md), [projective caps](PROJECTIVE_CAP_PROOF.md), [all-short constraints](ALL_SHORT_PROOF.md), [kernel bridge](CYCLIC_CLASS_BRIDGE.md), [collinear extraction](COLLINEAR_CYCLIC_CREDIT.md), [centered-gradient gate](MEANVALUE_CYCLIC_GATE.md), and [origin-ratio correction](ORIGIN_RATIO_PROOF.md). The exact ratio is 7/10; the latter proof supersedes ratio-one origin terms in generic component notes. Dimension-specific illustrations in those notes are not premises: target coefficients and constraints are checked anew.

The mod-2 framework extends [Zubrilina (2026)](https://arxiv.org/abs/2610.35830), combined with Kravatskiy's cyclic-kernel methodology and new lattice collinear/rank refinements. No worldwide priority claim is made. Prior unrestricted bounds, if included in comparison tables, are not assumed in this proof. No unverified lattice lower comparison is included.

Comparison: the [preceding unrestricted source](../../dim44/refinement_profile20261114/README.md) claims **<0.0000018326064498** (verified at commit c6f6ea00). It also bounds lattices. The decrease is approximately 2.6553%; this is comparison data, not a proof premise.

Endpoint convention: bins are right-open [E_j,E_{j+1}), and cumulative cutoffs are strict. Closed intervals are used only to verify conservative envelope inequalities. Thus A=0 when short_upper+long_upper≤4 is valid; vectors at the final squared norm 9/4 belong to the nonpositive tail. The j=−1 branch means a class with no short vector, not an empty class.
