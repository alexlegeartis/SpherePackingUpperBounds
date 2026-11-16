# 37D lattice-only bound

Every 37-dimensional lattice packing has density **<0.000018573668945**. This is not an unrestricted-packing claim.

Install Python 3.12+ and the pinned requirements. From the repository root run `python verifications/lattice/replay.py --dimension 37`; alternatively run `python run_all.py` inside a materialized package. No optimizer, network or prior density/classification assumption is needed. Assertions must remain enabled.

Certificate SHA256: `5b6ae7366e992d6af08549ca0c29ea4f4e6a45bb30f5101bf76029d41341661f`.

The complete centered-gradient angular gate is shared by the two independent radial/profile audits. Its analytic derivative derivation was independently reviewed; this is not a claim of two independent angular implementations. The runner replays the whole gate, then both independent radial/profile implementations, including all 240 bins, 0 AA cuts, 193 class branches, Fourier/global signs, weighted rank inequalities, exact envelopes and normalization. Both algorithms use Arb.

The new 37D angular polynomial has SHA e3b262f8cb8e9b4db4023ac33bf0b6c7e58afa1818a6c33b347f5d31253c7549 and its complete cyclic domain is rechecked.

See [class foundations](CLASS_PROFILES.md), [angular proof](ANGULAR_PROOF.md), [projective caps](PROJECTIVE_CAP_PROOF.md), [all-short constraints](ALL_SHORT_PROOF.md), [kernel bridge](CYCLIC_CLASS_BRIDGE.md), [collinear extraction](COLLINEAR_CYCLIC_CREDIT.md), [centered-gradient gate](MEANVALUE_CYCLIC_GATE.md), and [origin-ratio correction](ORIGIN_RATIO_PROOF.md). The exact ratio is 7/10; the latter proof supersedes ratio-one origin terms in generic component notes. Dimension-specific illustrations in those notes are not premises: target coefficients and constraints are checked anew.

The mod-2 framework extends [Zubrilina (2026)](https://arxiv.org/abs/2610.35830), combined with Kravatskiy's cyclic-kernel methodology and new lattice collinear/rank refinements. No worldwide priority claim is made. Prior unrestricted bounds, if included in comparison tables, are not assumed in this proof. No unverified lattice lower comparison is included.

Comparison uses the preceding lattice-only certificate in ../dim37/README.md, retained unchanged. The earlier unrestricted result in ../../dim37/refinement_profile20261114/README.md is separate, not a proof premise. No worldwide priority claim.

Radial bins are right-open; cumulative cutoffs are strict. Branch j=-1 has no short vector. The nonpositive tail includes9/4.
