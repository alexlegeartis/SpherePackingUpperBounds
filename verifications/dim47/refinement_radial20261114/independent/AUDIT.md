# Independent implementation audit

Fresh `audit_angular_cyclic.py` and `audit_radial.py` replays both pass for
dimension 47. `angular_result.json`, `cyclic_result.json` and
`radial_result.json` record their actual results. The frozen witness data are
copied under `inputs/` and hash-linked by the parent entrypoint.

The angular audit reconstructs the target-dimensional Gegenbauer basis by a
three-term recurrence, checks its constant coefficient by spherical moments,
uses rational Sturm chains, and proves the complete compact cyclic domain
with changed radicand enclosures. The scale and infinite-tail arguments are
given in [the proof](../PROOF.md).

The radial audit reconstructs the Fourier transform by recurrence, independently
of the primary coefficient-sum formula. It checks the full signed pair gap,
the exact union of all accepted intervals, the origin Gram identity, positive
normalization, and rational strict rounding. It shares FLINT/Arb root
isolation with the primary implementation. Independence here concerns code
and enclosures; it does not mean an independent arithmetic library or an
external human audit. No optimization solver is used.
