# Certified sphere-packing constructions beyond 36 dimensions

Research checkpoint, 13 November 2026. The principal result is an explicit
**65-dimensional even integral lattice** with

\[
\min L_{65}=12,\qquad \det L_{65}=10\cdot3^{31},\qquad
\boxed{\delta(L_{65})=\frac{3^{17}}{\sqrt{10}}
=40,837,705.2485403\ldots.}
\]

The crucial new computation is a complete, independently reviewed interval
enumeration proving that an explicit dual coset of Ne64 has minimum
`26/3`. All four shards finished, visiting **4,880,529,145 nodes** including
repeated upper-tree work. The saved 65D Gram and its exact minimum/determinant
proof are complete.

This checkpoint also gives nine exactly checked periodic packings in
dimensions **60–62 and 66–71**. These and the 65D lattice exceed the specific
comparison constructions audited here. **Global record status and publication
priority have not been established.** The public density tables are sparse
in these dimensions; no missing table entry is treated as evidence of novelty.
The existing record registry is not changed.

All densities below are **center densities**. The occupied-volume fraction is
`Delta_n = V_n delta_n`, with `V_n = pi^(n/2)/Gamma(n/2+1)`.
Determinants refer to Gram determinants, so lattice covolume is their square root.

## Main results

| Dimension | Construction | Certified center density | Comparison established here |
|---|---|---:|---|
| 60 | Five Ne64 antipode cosets | 773,597.146771 | Exceeds the audited propagated families and the lower bound printed in Vardy 1999 |
| 61 | Four Ne64 antipode cosets | 2,043,183.683386 | Same comparison scope |
| 62 | Three Ne64 antipode cosets | 5,188,971.717320 | Same comparison scope |
| **65** | **Explicit Ne64 deep-coset extension** | **40,837,705.248540** | **9.5445115% above the depth-8 extension established in this run** |
| 66 | Nine Gamma72 antipode cosets | 82,897,968.185303 | 6.978039% above the inherited-minimum bare-section ceiling |
| 67 | Six Gamma72 antipode cosets | 247,629,167.976493 | 30.459863% above that ceiling |
| 68 | Five Gamma72 antipode cosets | 812,141,667.469561 | 51.273174% above that ceiling |
| 69 | Four Gamma72 antipode cosets | 2,342,574,740.00186 | 54.268973% above that ceiling |
| 70 | Three Gamma72 antipode cosets | 7,120,856,007.29983 | 43.582984% above that ceiling |
| 71 | Two Gamma72 antipode cosets | 19,663,555,319.9372 | 14.456956% above that ceiling |

The Gamma72 ceiling applies to bare-section **lower bounds computed using
the inherited parent minimum 8**. It does not bound sections whose actual
minimum is higher, other parents, other antipode clusters, or arbitrary packings.
The Vardy comparisons concern his printed lower bounds, not an asserted exact
density for every Elkies section or a current SOTA table.

### Exact antipode parameters

Here `K` is the primitive normal Gram, `s` is the number of projected cosets,
and `beta` is their squared diameter in the normal projection.

| Dimension | Parent | det K | s | beta |
|---|---|---:|---:|---:|
| 60 | Ne64 | 125 | 5 | 2/5 |
| 61 | Ne64 | 1024/27 | 4 | 3/8 |
| 62 | Ne64 | 12 | 3 | 1/3 |
| 66 | Gamma72 | 12288 | 9 | 1/2 |
| 67 | Gamma72 | 13122 | 6 | 2/9 |
| 68 | Gamma72 | 2000 | 5 | 1/5 |
| 69 | Gamma72 | 300 | 4 | 1/5 |
| 70 | Gamma72 | 48 | 3 | 1/6 |
| 71 | Gamma72 | 10 | 2 | 1/10 |

For these certificates the exact formula is

\[
\delta_n=\frac{s}{\operatorname{covol}(\Lambda)\sqrt{\det K}}
\left(\frac{\mu-\beta}{4}\right)^{n/2}.
\]

The parent parameters are `(mu,covol)=(12,3^16)` for Ne64 and `(8,1)` for
Gamma72. In the Ne64 construction the normal lattice lies in the parent
dual. Exact saturation makes the evaluation map onto the integer label
lattice; the period covolume is `covol(parent)*sqrt(det K)`.

Exact checkers verify the published Gram, positive definiteness, determinants,
primitive normal embeddings, every cluster distance, and rational squared
densities. The Ne64 audit additionally constructs complete integer bases and
rational shifts of the actual 60–63D period lattices. The parent minima are
inherited from published theorems, explicitly recorded as such.

## The 65D construction and its proof

Let `L=Ne64`. Its published properties are minimum 12, determinant `3^32`,
and 3-modularity. An explicit vector `h` in `L*`, of order 3 modulo `L`, satisfies

\[
\min_{x\in L}\|h+x\|^2=\frac{26}{3}.
\]

Set `t=sqrt(10/3)` and define

\[
L_{65}=\{(x+mh,mt):x\in L,\ m\in\mathbb Z\}.
\]

The minimum proof covers all layer indices:

* For `m=0`, a nonzero vector has squared norm at least 12 by the parent theorem.
* For `m` not divisible by 3, the horizontal component belongs to `h+L` or
  `-h+L`, so the squared norm is at least `26/3+10/3=12`.
* For nonzero `m` divisible by 3, the vertical squared norm is at least
  `9*(10/3)=30`.

The saved generator itself has norm 12. Thus the minimum is exactly 12.
The Gram is even integral, and its Schur complement over the parent Gram
is `10/3`. Consequently

\[
\det L_{65}=3^{32}\frac{10}{3}=10\cdot3^{31},\qquad
\delta(L_{65})=\frac{(\sqrt{12}/2)^{65}}{\sqrt{10\cdot3^{31}}}
=\frac{3^{17}}{\sqrt{10}}.
\]

### What was exhaustively certified

In an exactly checked reduced basis with Gram `Q`, the numerator of every
vector in the selected dual coset has coordinates `z=c (mod 3)` and norm
`zQz^t=6 (mod 18)`. We exhaustively excluded `zQz^t<=60` in that residue
class. The next allowed value is 78; a saved exact norm-78 witness proves
the coset minimum is exactly `78/9=26/3`.

The independent enumerator uses an exact rational LDL decomposition enclosed
in dyadic intervals, exact 128-bit center-endpoint sums, outward-rounded
binary64 bounds for pruning, and exact integer leaf checks. Four deterministic
shards cover every prefix. All finished with no violating leaf. Their input,
source, and run hashes are retained. A second full traversal with an exact
global rounding-error audit is included as a supplementary check.

The primary proof, the arithmetic assumptions, and the complete reproduction
command are in [ne64_hole_certificate/README.md](ne64_hole_certificate/README.md).
The saved matrix is [higher_parent/ne65_candidate_gram.txt](higher_parent/ne65_candidate_gram.txt);
the historical word `candidate` in the filename is not its current proof status.
[higher_parent/ne65_verified.json](higher_parent/ne65_verified.json) records
`verified_minimum_12_lattice`.

### Exact comparison in 65D

The theta series and the discriminant-class count give a separate existential
Ne64 coset of minimum at least 8. That yields the baseline

\[
\delta_{65}\ge\frac{3^{33/2}}2=37,279,553.9356211\ldots.
\]

Our explicit coset improves this by the exact factor

\[
\sqrt{\frac65}=1.095445115\ldots.
\]

This is a comparison between two rigorously established constructions. It is
not a statement that the depth-8 construction was the previously published
world record.

## Benchmarks that disqualified apparent advances

The initial Ne64 sections in 57–59D and antipode in 63D are retained as valid
constructions and controls, but are excluded from the main candidate list.

| Dimension | Initial construction | Stronger comparison found |
|---|---:|---|
| 57 | 57,530.175827 | A torus-layer consequence of the published Schuett56 lattice exceeds 151,706 |
| 58 | 115,060.351653 | A primitive Gamma72 section gives `2^18/sqrt(3)=151,348.908966...` |
| 59 | 244,079.864699 | A primitive Gamma72 section gives `2^18=262,144` |
| 63 | 13,105,940.958282 | Vardy 1999 already reports an Elkies lower bound with `log2(delta)=23.92...` |

The 58D/59D comparisons are certified using the explicit scaled-Leech block
of Gamma72 and primitive laminated normal embeddings. The audit also checks
products, torus colorings, and averaged ball-window projections from the
published 56D, 64D and 72D parents. All families and their limitations appear
in [audits/PROPAGATED_BENCHMARKS.md](audits/PROPAGATED_BENCHMARKS.md).

The current public Cohn table covers dimensions 1–48 and 56, 64, 72. The
Nebe–Sloane high-dimensional density table is also sporadic. Targeted searches
did not establish an exhaustive current comparison for 60–62 or 65–71.
See [code_routes/BENCHMARK_AUDIT.md](code_routes/BENCHMARK_AUDIT.md).

## Exclusions and additional mathematical results

These results constrain further searches; they are not density records.

* **71D rank-one optimum within the chosen formula.** Among every rank-one
  Gamma72 antipode bound using parent minimum 8, the two-coset construction
  with normal norm 10 is optimal. This includes every number of labels and
  every allowed even normal norm, not just the enumerated examples.
* **Unimodular-neighbor obstruction.** In any even unimodular 72D lattice
  containing a vector of norm at most 6, every nonzero lattice vector of norm
  at most 10 has a vector of norm at most 6 in its orthogonal section. A degree-two
  harmonic theta identity proves this, ruling out the proposed short-normal
  neighbor route to minimum-8 sections.
* **71D determinant corollary.** Every even integral positive-definite 71D
  lattice of minimum at least 8 has determinant at least 8. The proof combines
  that obstruction with exact small discriminant-form classification and gluing.
* **33D layer-family exclusion.** For each of seven catalogue32D parents,
  all 65,536 dual quotient classes have explicit short representatives. A
  finite quadratic-space argument proves that arbitrary heights, varying
  overall minimum, and arbitrary sequences of dual-offset layers cannot
  beat `3^(33/2)/2^25` in 33D. The proof is restricted to this stated family.
* **Rejected order-three shortcut.** The known Ne64 automorphism subgroup
  acts faithfully on its discriminant group. The fixed-point-free order-three
  elements in its SL2(5) factor fix only isotropic classes, so they do not prove the desired
  nonisotropic depth-26/3 hole by simple orbit counting. The successful 65D
  result uses the explicit exhaustive coset computation instead.

Additional exploratory Gram relaxations in dimensions 33 and 34 remained
far below their benchmarks; their conditional floating outputs are retained
only as search history.

## Reproduction and directory guide

Run the compact algebraic audit, including the stored interval completion
records and all finite covering witnesses:

```bash
OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/verify_checkpoint.py
```

To repeat the complete 65D supporting enumeration:

```bash
OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/ne64_hole_certificate/run.py
```

The compact checks require Python, NumPy, SymPy, and mpmath. The full interval
replay also requires g++. Optional discovery programs use fpylll/SciPy.

| Directory | Contents |
|---|---|
| `ne64_hole_certificate/` | Primary exhaustive hole certificate and replay |
| `higher_parent/` | Explicit 65D lattice; Gamma72 antipodes; theta count; supplementary rounding audit |
| `dual64/` | Ne64 normal embeddings and finite antipode clusters |
| `audits/` | Independent period bases, exact checks, and propagated comparisons |
| `code_routes/` | Harmonic-theta obstruction, rank-one proof, discriminant actions, and priority audit |
| `transfer_neighbors/` | Seven complete32D quotient witnesses and33D layer obstruction |
| `low_contact_seeds/` | Unsuccessful, explicitly conditional Gram-relaxation history |

## Primary sources

1. [Ne64 catalogue and its minimum theorem references](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ64b.html).
2. [Gamma72 catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html) and Nebe's [minimum-8 construction paper](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/papers/Ext72.pdf).
3. Conway and Sloane, [The Antipode Construction for Sphere Packings](https://neilsloane.com/doc/anti.pdf).
4. Chen, Hu, Li, Wang and Wu, [New Sphere Packings from the Antipode Construction](https://arxiv.org/abs/2505.02394), 2025; its stated dimensions are 19, 20, 21, 23, 44, 45, 47.
5. Vardy, [Density Doubling, Double-Circulants, and New Sphere Packings](https://citeseerx.ist.psu.edu/document?doi=c9085ab03f7ab8f5154d3c38fe2be7ab6316fac4&repid=rep1&type=pdf), Trans. AMS351(1999),271–283, section4.
6. [Cohn's current density table](https://cohn.mit.edu/sphere-packing/) and [Nebe–Sloane density table](https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html).

Specific modular-form, theta-series, and discriminant-form references are
given beside the corresponding proof in `code_routes/`.
