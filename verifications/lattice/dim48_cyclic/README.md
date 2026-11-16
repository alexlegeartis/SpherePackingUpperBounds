# 48D lattice-only bound

Every 48-dimensional lattice packing has density **<0.000000358368032**. This is not an unrestricted-packing claim.

Install Python 3.12+ and the pinned requirements. From the repository root run `python verifications/lattice/replay.py --dimension 48`; alternatively run `python run_all.py` inside a materialized package. No optimizer, network or prior density/classification assumption is needed. Assertions must remain enabled.

Certificate SHA256: `a03938f8d309fcc2c237bd0df6342f45c44fced0765a1b453a66f85f5f6bad6f`.

The complete centered-gradient angular gate is shared by the two independent radial/profile audits. Its analytic derivative derivation was independently reviewed; this is not a claim of two independent angular implementations. The runner replays the whole gate, then both independent radial/profile implementations, including all 240 bins, 0 AA cuts, 193 class branches, Fourier/global signs, weighted rank inequalities, exact envelopes and normalization. Both algorithms use Arb.

The new 48D angular polynomial has SHA c89a694a264a3ebdcff9690a12aba03bca7e6be14c768d5465e10e146c914700 and its complete cyclic domain is rechecked.

See [class foundations](CLASS_PROFILES.md), [angular proof](ANGULAR_PROOF.md), [projective caps](PROJECTIVE_CAP_PROOF.md), [all-short constraints](ALL_SHORT_PROOF.md), [kernel bridge](CYCLIC_CLASS_BRIDGE.md), [collinear extraction](COLLINEAR_CYCLIC_CREDIT.md), [centered-gradient gate](MEANVALUE_CYCLIC_GATE.md), and [origin-ratio correction](ORIGIN_RATIO_PROOF.md). The exact ratio is 7/10; the latter proof supersedes ratio-one origin terms in generic component notes. Dimension-specific illustrations in those notes are not premises: target coefficients and constraints are checked anew.

The mod-2 framework extends [Zubrilina (2026)](https://arxiv.org/abs/2610.35830), combined with Kravatskiy's cyclic-kernel methodology and new lattice collinear/rank refinements. No worldwide priority claim is made. Prior unrestricted bounds, if included in comparison tables, are not assumed in this proof. No unverified lattice lower comparison is included.

Comparison: preceding lattice table upper bound 0.000000408175; preceding unrestricted upper bound 0.00000040843570355 (SpherePackings main 7e899e58). These are distinct comparison data, not proof premises. The preceding lattice-only proof is retained at ../dim48/README.md; unrestricted source: ../../dim48/refinement_larger20261115/README.md. The existing lattice-table lower comparison 2.31782953e-8 is retained.

Endpoint convention: radial bins are right-open [E_j,E_{j+1}) and cumulative cutoffs are strict. Closed intervals verify conservative envelopes. Vectors at9/4 belong to the nonpositive tail. Branch j=-1 denotes a class with no short vector, not a literally empty class.

Its own dimension48 angular gate is replayed completely; no dimension-transfer premise is used.
