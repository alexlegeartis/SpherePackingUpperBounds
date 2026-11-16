# 43D lattice-only bound

Every 43-dimensional lattice packing has density **<0.000002566728865**. This is not an unrestricted-packing claim.

Install Python 3.12+ and the pinned requirements. From the repository root run `python verifications/lattice/replay.py --dimension 43`; alternatively run `python run_all.py` inside a materialized package. No optimizer, network or prior density/classification assumption is needed. Assertions must remain enabled.

Certificate SHA256: `6c439b6ee1d2afae624f21587c8ad79cc715985b741b4c7b9c1234e67cb86239`.

The complete centered-gradient angular gate is shared by the two independent radial/profile audits. Its analytic derivative derivation was independently reviewed; this is not a claim of two independent angular implementations. The runner replays the whole gate, then both independent radial/profile implementations, including all 240 bins, 2 AA cuts, 193 class branches, Fourier/global signs, weighted rank inequalities, exact envelopes and normalization. Both algorithms use Arb.

The new43D angular polynomial has SHA12967ed055a9b6890a23771e9e737927bcd228bb7f7fc1a1530026005f8aa62e and its complete cyclic domain is rechecked.

See [class foundations](CLASS_PROFILES.md), [angular proof](ANGULAR_PROOF.md), [projective caps](PROJECTIVE_CAP_PROOF.md), [all-short constraints](ALL_SHORT_PROOF.md), [kernel bridge](CYCLIC_CLASS_BRIDGE.md), [collinear extraction](COLLINEAR_CYCLIC_CREDIT.md), [centered-gradient gate](MEANVALUE_CYCLIC_GATE.md), and [origin-ratio correction](ORIGIN_RATIO_PROOF.md). The exact ratio is 7/10; the latter proof supersedes ratio-one origin terms in generic component notes. Dimension-specific illustrations in those notes are not premises: target coefficients and constraints are checked anew.

The mod-2 framework extends [Zubrilina (2026)](https://arxiv.org/abs/2610.35830), combined with Kravatskiy's cyclic-kernel methodology and new lattice collinear/rank refinements. No worldwide priority claim is made. Prior unrestricted bounds, if included in comparison tables, are not assumed in this proof. No unverified lattice lower comparison is included.

Comparison: the [preceding unrestricted source](../../dim43/refinement_profile20261114/README.md) claims **<0.0000026096189902** (verified at commitc6f6ea00). It also bounds lattices. The decrease is approximately1.6435%; this is comparison data, not a proof premise.

Endpoint convention: radial bins are right-open [E_j,E_{j+1}) and cumulative cutoffs are strict. Closed intervals only verify conservative envelopes. The equality test A=0 at short_upper+long_upper≤4 is valid. Vectors at9/4 belong to the nonpositive tail. Branch j=−1 denotes a class with no short vector, not a literally empty class.
