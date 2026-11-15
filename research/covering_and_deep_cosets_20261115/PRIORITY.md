# Priority verdicts for the five additions from e9d03c1

**Verdict:** the previous five values are superseded by stronger verified
bounds in this continuation. They should not be described as current
worldwide records. The following verdicts concern the actual numbers,
not a renamed status category.

| Dimension | Earlier center density | Stronger verified replacement | Current verdict on earlier value |
| --- | ---: | ---: | --- |
| 49 | 13,039.647236681… | 14,578.768811918… | Superseded by an explicit deeper-coset lattice |
| 50 | 11,663.015049534… | 12,625.584147021… | Superseded by an explicit fitted-metric lattice |
| 51 | 11,064.507582405… | 11,794.804713454… | Superseded by an explicit fitted-metric lattice |
| 52 | 12,120.560781140… | 12,370.495545910… | Superseded by an explicit fitted-metric lattice |
| 81 | 996,397,626,643.825… | 1,066,248,855,605.296… | Superseded by the stronger covering-based consequence of MW80 |

Every replacement is checked against the pinned rational squared
density in [previous_bounds.json](previous_bounds.json). In 49–52D,
the explicit parent, residue rows and full minimum proof are available.
In 81D, the full existence argument is supplied, without explicit shifts.

## Bibliographic evidence actually recovered

* **Cohn and Nebe–Sloane record tables.** Their coverage is sparse beyond
  48D. Cohn's table includes 56D, 64D and 72D, while the catalogue and
  the SPLAG third-edition preface add selected parent dimensions. None
  of those inspected tables identifies a previous best result in all
  five target dimensions. Absence is not a record confirmation.
* **Classical P48 extension methodology.** BCS1982 supplies the P48q
  parent and its norm-doubler. SPLAG's third-edition preface, p.xix,
  describes Elkies's extension of Q32 using half-lattice holes of norm
  4. This is relevant prior construction methodology, not a reported
  numerical 49–52D comparator. We do not claim to originate the
  half-coset/lamination principle.
* **Shioda's MW80.** Schütt–Shioda's Theorem 13.35 at p=41 gives the
  exact parent density `7^40/(2^40*41^6)`. It does not print the previous
  81D value as a worldwide record. The newly recovered covering theorem
  yields a strictly stronger 81D consequence of that same parent.
* **Coxeter–Few–Rogers (1959).** The publisher's original abstract
  confirms the universal covering lower bound. Fejes Tóth's 2022
  author survey, §2.2, gives the precise regular-simplex statement.
  The finite rational estimate and its application to torus layers
  are derived in [COVERING_PROOF.md](COVERING_PROOF.md). The full
  original 1959 article was not downloaded; its theorem is inherited.
* **Other examined construction papers.** Chen2013, Shimada2012,
  Leterrier2022, Nigam2022 and the earlier documented Mordell–Weil and
  antipode sources do not supply a stronger, reported numerical bound
  in these five dimensions in the inspected statements and tables.
  Leterrier's 54D guarantee materially strengthens the nearby 53D
  comparison; it is used to reject weaker 53D proposals.

Queries included dimension names written out and numerically, P49–P52,
P48/P48q holes and covering radii, MW54/MW80 and their sections,
laminated extensions, and current record tables. The search did not
recover a trusted five-dimensional previous-record list that would
justify a worldwide SOTA claim for the replacements. We therefore
confirm their construction proofs and specified external comparisons,
and leave their bibliographic priority unconfirmed. This is an actual
missing source attribution, not a numerical-certificate failure.

## The separately confirmed 96D comparison

BCS1982 **explicitly reports** its P48q-based 96D construction in
Table 1, p.177. Equation (15) reconstructs the exact center density
`3^48/2^24`. The stronger current bound has center density
`10,391,088,344,672,478.938…`, giving the verified gain
**118.554979900061…% over that reported construction**.
The exact rational-square comparison is in
[covering_results.json](covering_results.json).

That source attribution and gain are confirmed. It is not an assertion
that a search can prove the absence of all other unpublished or
unlocated constructions. No statement of global optimality is made.
