# Density benchmark and priority audit,13 November2026

## Scope and conclusion

The current primary public tables are sparse above dimension48. The
search below did not establish a complete best-known-density benchmark
for every dimension57–71. Therefore this directory makes no global SOTA
or publication-priority assertion for the new constructions.

The substantive correction found here is the newer Schuett56D input.
Simple propagation from it already invalidates a proposed57D record
comparison, and a torus-layer propagation invalidates the initial58D
comparison. The other research agent subsequently found stronger58D and
59D comparisons from Gamma72's Leech sections. This illustrates why an
old table row alone is not an adequate baseline.

## Primary sources checked

|Source|What it establishes|What it does not establish|
|---|---|---|
|Henry Cohn, current sphere-packing table, <https://cohn.mit.edu/sphere-packing/>|Exact current rows through48 and at56,64,72; row56 is5^28/(2^28*29^4)|Rows for57–63 or65–71|
|Nebe–Sloane catalogue density table, <https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html>|Same principal high-dimensional anchors;56D names Schuett, Mordell-Weil Lattices2019,p400|An exhaustive list of derived sections/laminations between anchors|
|Schuett–Shioda, *Ranks and Sphere Packings*,2019, <https://doi.org/10.1007/978-981-32-9301-4_13>|Primary chapter identity and publication date; linked from the56D table reference|The chapter's full text was not available in the inspected publisher preview|
|Nebe, *An even unimodular72-dimensional lattice of minimum8*, <https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/papers/Ext72.pdf>|Gamma72's construction, minimum8, and explicit scaled-Leech sublattice structure|A complete density table for its lower-dimensional sections or antipodes|
|Conway–Sloane, *The antipode construction for sphere packings*, <https://neilsloane.com/doc/anti.pdf>|The general antipode construction and its published lower-dimensional examples|No64D/72D-derived numerical benchmark appears in the inspected paper|
|Recent authors' reproducible repository, <https://github.com/wcgbg/sphere-packing>|Current work on extremal48D sections and antipodes, including38–43 density claims|No benchmark for58–71 was identified there|

The catalogue HTML downloaded in this audit is retained as
`density_catalogue.html`. Its header still says February2012 while its
56D row cites2019; that header must not be treated as a reliable date for
all entries.

## Concrete propagated comparison bounds

Write `d56=5^28/(2^28*29^4)=196213.21381197087...`.

* Direct product with a one-dimensional unit-radius packing gives
  `delta57 >= d56/2 =98106.6069059854...`. This alone exceeds the initial
  Ne64 bare-section candidate57530.1758....
* The exact torus-layer certificate in `TORUS_LAYER_BENCHMARKS.md` gives
  `delta57 >= d56 sqrt(1000/1673) =151698.1486935438...`.
* Its triangular three-color version gives
  `delta58 >= 200d56/(173 sqrt(3)) =130963.8749282115...`, exceeding the
  initial Ne64 bare-section candidate115060.3516....
* The coordinating team's exact Gamma72/Leech-section audit gives still
  stronger comparisons `delta58 >=2^18/sqrt(3)` and `delta59 >=2^18`.
  Those proofs and matrices are in the sibling audit/transfer directories,
  not duplicated here.

The torus bounds are derived consequences of a published parent. Their
specific numerical consequences were not located as separately published
rows, so they should not be falsely attributed to Schuett or Cohn.

## Priority limits for60–63 and65–71

Targeted searches for Ne64 and Gamma72 together with sections, projections,
antipodes, and the individual dimensions found the sources above but no
comprehensive current table in these ranges. Searches by spelled-out
dimension and Hermite-constant terminology did not resolve the gap.

The older Vardy1999 section bounds are useful historical comparisons, as
is the norm8 Gamma72 hyperplane-section value. They are not a proof that
all later constructions or routine consequences have been exhausted.
The current candidates should therefore be described as certified lower
bounds, with precisely stated comparison families. A global record claim
requires a stronger external priority audit than this search supplied.
