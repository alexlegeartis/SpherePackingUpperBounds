# 48D lattice-only bound

Every 48-dimensional lattice packing has density **<0.000000354687644**. This is not an unrestricted-packing claim.

Install Python 3.12+ and the pinned requirements. From the repository root run `python verifications/lattice/replay.py --dimension 48`; alternatively run `python run_all.py` inside a materialized package. No optimizer, network or prior density/classification assumption is needed. Assertions must remain enabled.

Certificate SHA256: `a40d47f8297b6475f8173ae8ba3f1c84356ecc36c43c18094b3cf31dc3c85df1`.

The complete centered-gradient angular gate is shared by the two independent radial/profile audits. Its analytic derivative derivation was independently reviewed; this is not a claim of two independent angular implementations. The runner replays the whole gate, then both independent radial/profile implementations, including all 120 bins, 0 AA cuts, 97 class branches, Fourier/global signs, weighted rank inequalities, exact envelopes and normalization. Both algorithms use Arb.

The new 48D angular polynomial has SHA c89a694a264a3ebdcff9690a12aba03bca7e6be14c768d5465e10e146c914700 and its complete cyclic domain is rechecked.

See [class foundations](CLASS_PROFILES.md), [angular proof](ANGULAR_PROOF.md), [projective caps](PROJECTIVE_CAP_PROOF.md), [all-short constraints](ALL_SHORT_PROOF.md), [kernel bridge](CYCLIC_CLASS_BRIDGE.md), [collinear extraction](COLLINEAR_CYCLIC_CREDIT.md), [centered-gradient gate](MEANVALUE_CYCLIC_GATE.md), and [origin-ratio correction](ORIGIN_RATIO_PROOF.md). The exact ratio is 7/10; the latter proof supersedes ratio-one origin terms in generic component notes. Dimension-specific illustrations in those notes are not premises: target coefficients and constraints are checked anew.

The mod-2 framework extends [Zubrilina (2026)](https://arxiv.org/abs/2610.35830), combined with Kravatskiy's cyclic-kernel methodology and new lattice collinear/rank refinements. No worldwide priority claim is made. Prior unrestricted bounds, if included in comparison tables, are not assumed in this proof. No unverified lattice lower comparison is included.

Comparison uses the preceding lattice-only certificate in ../dim48_cyclic/README.md, retained unchanged. The earlier unrestricted result in ../../dim48/refinement_larger20261115/README.md is separate, not a proof premise. No worldwide priority claim.

Radial bins are right-open; cumulative cutoffs are strict. Branch j=-1 has no short vector. The nonpositive tail includes9/4.
