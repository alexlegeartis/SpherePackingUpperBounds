# Stronger bounds from covering geometry and deeper compatible cosets

**Current continuation:** [adaptive layers](../adaptive_layers_20261115/) strengthens13 bounds, including50D–53D,77D and85D by9.99–18.53%. The checkpoint numbers below are superseded where indicated; the [current comparison table](PRIORITY.md) uses the new results.

**All five numbers from the previous update are now superseded.** The
earlier 49–52D and 81D values should not be presented as current records:
the verified constructions below are strictly stronger. The additional
searches also improve **57D, 82D and 93–96D**.

All numbers are **center densities**, not occupied fractions. Displayed
decimals are approximate; the reports retain exact rational squares. Changes in
the last column compare with the pinned repository checkpoint
`e9d03c1`; they are internal progress, not percentages against a
previously reported worldwide record.

| Dimension | Current bound | Method | Increase over checkpoint |
| --- | ---: | --- | ---: |
| **49** | **14,578.7688119179…** | Explicit P48q extension, one deeper coset | **11.803399%** |
| **50** | **12,625.5841470212…** | Explicit compatible cosets and fitted vertical Gram | **8.253175%** |
| **51** | **11,794.8047134541…** | Same, with exact vertical determinant 55/8 | **6.600358%** |
| **52** | **12,370.4955459097…** | Same, with exact vertical determinant 75/8 | **2.062073%** |
| **57** | **163,056.783695743…** | Covering bound + MW56, two colors | **7.482034%** |
| **81** | **1,066,248,855,605.2958…** | Covering bound + MW80, two colors | **7.010377%** |
| **82** | **1,040,071,302,182.4955…** | Covering bound + MW80, three colors | **2.400937%** |
| **93** | **8.742047493270432e15** | Covering bound + MW92, two colors | **6.837494%** |
| **94** | **8.583004659539084e15** | Covering bound + MW92, four colors | **11.435576%** |
| **95** | **8.865551293806628e15** | Covering bound + MW92, six colors | **19.636193%** |
| **96** | **1.039108834467248e16** | Covering bound + MW92, nine colors | **26.156816%** |

The new **96D** result is **118.554980% above the explicitly reported
BCS1982 construction** `delta=3^48/2^24`. Its occupied fraction is
approximately **6.108662600212e-22**. This is the verified external gain;
it supersedes the earlier 73.240722% comparison to the same source.

There are also sharper **derived section controls in 53D, 55D and
89D–91D**. A finite estimate of Rogers's packing simplex bound improves
the previous Blichfeldt Hermite estimates. For example the new 89D
guarantee is approximately **7.144650443536e13**, and the 53D guarantee
is approximately **17,624.661756**. These are conservative consequences
of published parents, not new explicit normals or reported records.
[Derivation and all five values](SECTIONS_PROOF.md) · [exact audit](verify_sections.py).

The [remaining-dimensions review](REMAINING_DIMENSIONS.md) extends the
best-known claims to33 construction dimensions and corrects the 55D
literature entry to Leterrier2024's **116,783.149646…**, above our old
55D control. The explicit53D lattice now replaces that section estimate;89D–91D remain classical derived controls.

## What the priority investigation confirms

**The bounds in 49D–53D, 57D, 59D–62D, 65D–71D, 73D–79D, 81D–85D and 93D–96D are best-known
lower bounds in our literature review dated 15 November 2026.** No
stronger published lower bound was located. [Sources and comparison table](PRIORITY.md)
identify the work and theorem behind each numerical comparison.

We compare proved lower bounds: B>A improves a paper's guarantee A.
Neither a previous-record table nor the unknown actual minimum of an
older lattice is required. The earlier five repository values are
superseded by the current constructions. The percentages in the first
table above track that internal progress; the external gains are in
[PRIORITY.md](PRIORITY.md) and [literature_comparisons.json](literature_comparisons.json).

The [full table for every dimension 48–96](../range48_96_20261115/TABLE.md)
uses the same dated review and distinguishes reported numbers from
derived consequences. Exact comparison arithmetic can be replayed with
`python research/covering_and_deep_cosets_20261115/verify_literature_comparisons.py`.

## Two ideas that produced gains

**Use a covering obstruction when selecting offsets.** If a few balls
covered the horizontal torus, their periodic lifts would cover Euclidean
space. Coxeter–Few–Rogers's classical covering-density theorem rules
that out at a substantially larger radius than the old volume-only
argument. A finite Jensen estimate gives the exact usable bound
`C_n=(n+1)*(n/(n+3))^(n/2)`. Applying it to all 200 previously tested
parent/color combinations produces the seven gains above.
[Complete derivation](COVERING_PROOF.md) · [exact checker](verify_covering.py) ·
[density squares and exact gain intervals](covering_results.json).

**Fit the vertical metric to different horizontal depths.** One P48q
half-coset is certified to depth 4, while the other classes only need
depth 3 or 7/2. A uniform root-lattice scaling wastes that distinction.
The new rational metrics exploit it. The full release replay checks
all fifteen required classes, with **638,998,811 enumeration nodes**,
and an exact-witness negative control. The small vertical inequalities
and full Gram determinants are checked exactly.
[Explicit proof and Gram formulas](DEEPER_PROOF.md) ·
[integer residues and density squares](deeper_results.json) ·
[verifier](verify_deeper.py).

The 49–52D shifts are explicit. The covering-based shifts in the other
seven dimensions remain existential. These methods use classical
construction ingredients; their present numerical consequences and
certificates are computed here.

## Other dimensions: tested ideas and next targets

* **53–55D, weighted layers.** A weighted Voronoi search found the
  successful 52D metric. Extending it with eight tested fifth-layer
  quadratic patterns produced no score above the classical 53D
  section control. Those numerical candidates were not promoted to
  packing certificates. A broader search over four-dimensional metrics
  remains a concrete route to 53D, rather than repeating the same root
  lattice scale.
* **Deeper compatible subspaces.** A totally singular binary subspace
  of rank five, with every nonzero half-coset at depth 4, would give an
  attractive explicit 53D extension. The first deeper class was found,
  but six tested rank-two extension candidates all failed compatibility:
  at least one of the added class or its sum had a norm-12 witness.
  This is a bounded unsuccessful search, not a nonexistence proof.
  A next search can use verified parent automorphisms to preserve one
  deep class and filter short witnesses before complete enumeration.
* **58–62D and 89–91D, structured normals.** The covering scan does not
  beat the existing section/antipode controls there. The new Rogers
  estimates already sharpen the conservative 89–91D controls. Actual short dual
  normals in MW56 or MW92 could improve those section bounds. The
  missing work is to construct and saturate those normals, then certify
  the period determinant and cluster diameter; an abstract Gram or a
  Hermite estimate is not an embedding certificate.

[Discovery budgets and rejected candidates](discovery_summary.json) are
retained separately from positive proof data. No construction above
96D is used as an input.

## Replay

Standard-library covering arithmetic:

```sh
python research/covering_and_deep_cosets_20261115/verify_covering.py
python research/covering_and_deep_cosets_20261115/verify_sections.py
```

Exact setup and saved completed enumeration reports (NumPy, SymPy and
python-flint), followed by a fresh full replay with GNU C++:

```sh
python research/covering_and_deep_cosets_20261115/verify_deeper.py
python research/covering_and_deep_cosets_20261115/verify_deeper.py --output /tmp/deeper_p48_replay
```

The output directory must be new. The full enumeration has the stated
outward-rounded IEEE-754/compiler assumptions; the quick audit does not
rerun it. Parent-minimum and covering theorems are stated external
mathematical inputs. Sources, hashes, query scope and source-access
limitations are recorded in [sources.json](sources.json) and the linked
earlier surveys.
