# Certified packings from Gamma72 and a deep Ne64 coset

This checkpoint gives explicit, certified sphere-packing lower bounds in dimensions **65–71**. Its strongest 65-dimensional construction is an explicit extension of Ne64 constructed here through a dual coset of minimum `26/3`; a complete outward-interval enumeration certifies that minimum. The 66–71 constructions use finite dual clusters in low-rank normal sublattices of Nebe's `Gamma72`. Published minimum theorems for both parent lattices are inherited, rather than re-enumerated.

**World-record status has not been established.** The standard Cohn and Nebe–Sloane density tables do not tabulate dimensions 65–71. The Gamma72 comparisons below are against a rigorous ceiling for the *bare-section bounds that use the inherited parent minimum 8*, not against an independently established list of all previous packings. A section with higher actual minimum is outside that comparison.

## Explicit 65-dimensional lattice: completed certificate

The lattice in `ne65_candidate.npz` and `ne65_candidate_gram.txt` is now **certified to have minimum 12**, with an even integral Gram and exact determinant

\[
10\cdot3^{31}=6,176,733,962,839,470.
\]

The historical `candidate` filenames are retained for continuity; their proof status is recorded in `ne65_verified.json`. The resulting center density is

\[
\boxed{\delta_{65}=\frac{3^{17}}{\sqrt{10}}
=40,837,705.2485403\ldots.}
\]

It improves on the depth-8 Ne64 extension proved below by the factor `sqrt(6/5)`, or **9.5445115%**, and improves on this checkpoint's Gamma72 65D construction. These are explicit, stated comparisons; global priority remains unestablished.

The construction uses an exact dual vector `h` with squared norm `26/3` and the generator `(h,sqrt(10/3))` above Ne64. Its third multiple has horizontal component in Ne64. All nonzero layer classes have squared norm at least `26/3+10/3=12`; nonzero layers divisible by three have vertical squared norm at least 30. The saved 65D Gram has final diagonal 12 and small integer off-diagonal entries.

The independent proof in `../ne64_hole_certificate/` excludes all integer vectors `z` of norm at most 60 in the prescribed residue class modulo 3. The residue forces numerator norm `6 mod 18`, so the next possible norm is 78, which becomes `26/3` after division by 9. An exact norm-78 witness proves equality.

All **four outward-interval shards completed**, with **4,880,529,145 total nodes**, including repeated upper-tree nodes, and zero exact leaves. Their identical shared-prefix count was **42,255,279**, and the assigned counts cover that set exactly. The complete proof summary and input/source hashes are in `../ne64_hole_certificate/summary.json`. This exhaustive computation is the primary certificate. Capped exploratory searches are not used as proof.

`verify_ne65.py` checks the original catalogue Gram, the basis change, the exact coset witness, all residue conditions, the full 65D Gram, its determinant and Schur complement, the layer inequalities, and the independent certificate's input hashes. It reports `verified_minimum_12_lattice` after consuming the completed summary.

```sh
python research/novel_routes_20261113/higher_parent/verify_ne65.py
```

A separate logged floating traversal visited 4,567,808,611 nodes and exhausted the same coset. The exact rational audit in `audit_ne64_rounding.py` gives a pruning-error bound below `1.016e-8` against a radius inflation of approximately `1e-6`. `NE64_CVP_CERTIFICATE.md` explains this independent cross-check; the outward-interval certificate remains primary.

## Best certified Gamma72 bounds

For a rank-`k` primitive normal Gram `K`, a finite set of dual integer labels `S`, squared diameter `beta`, and `n=72-k`, the lower bound is

\[
\delta_n=\frac{|S|}{\sqrt{\det K}}\left(\frac{8-\beta}{4}\right)^{n/2}.
\]

| Dimension | Normal determinant | Cosets | Squared diameter | Center density | Gain over inherited-minimum section ceiling |
|---|---:|---:|---:|---:|---:|
| 71 | 10 | 2 | 1/10 | 19,663,555,319.9372 | 14.456956% |
| 70 | 48 | 3 | 1/6 | 7,120,856,007.29983 | 43.582984% |
| 69 | 300 | 4 | 1/5 | 2,342,574,740.00186 | 54.268973% |
| 68 | 2000 | 5 | 1/5 | 812,141,667.469561 | 51.273174% |
| 67 | 13122 | 6 | 2/9 | 247,629,167.976493 | 30.459863% |
| 66 | 12288 | 9 | 1/2 | 82,897,968.1853032 | 6.978039% |
| 65 | 32768 | 9 | 1/2 | 37,073,098.4117906 | 10.486443% |

The 67–71 certificates are `gamma72_simplex_rank1.npz` through `gamma72_simplex_rank5.npz`. The 66- and 65-dimensional certificates are `gamma72_root_cluster_rank6.npz` and `gamma72_root_cluster_rank7.npz`. Each contains the exact ambient Gram `G`, integer coefficient rows `B` of the normal basis, its Gram `K=B G B^T`, and the dual integer labels `S`.

The 65-dimensional Gamma72 value is superseded by the certified Ne64 extension above; it is retained as a construction and control.

The 64-dimensional control, from a 9-point cluster in `(4E8)^*`, gives only `19,144,465.6988231`, below the published Ne64 value `3^16=43,046,721`. It is **not** an improvement in dimension 64.

The seven-coset rank-6 simplex also gives a valid 66-dimensional bound `77,604,025.3019868`, with determinant `89383` and diameter `180/791`. It is superseded here by the nine-coset `4E6` cluster.

## Why the certificates define packings

Let `Gamma` be the published even unimodular parent, let `U` be the real span of the rows `B`, and let `V=U^perp`. Exact Smith invariant factors equal to 1 show that `K=Gamma intersect U` is primitive. Unimodularity gives

\[
\pi_U(\Gamma)=K^*,\qquad \det(\Gamma\cap V)=\det K.
\]

For each label `s` in `S`, take the projection to `V` of all lattice points whose projection to `U` has dual coordinates `s`. This is a translate of the period lattice `Gamma intersect V`. The checker verifies

\[
(s-t)K^{-1}(s-t)^T\le\beta<8
\]

for every pair of labels. Distinct labels give distinct translates: otherwise their ambient difference would be a nonzero vector of `K` with squared norm below 8. Every difference between distinct projected centers therefore has squared length at least `8-beta`, by the orthogonal decomposition and the parent minimum theorem. The period covolume is `sqrt(det K)`, proving the stated lower bound.

The exact checker additionally re-parses the archived catalogue source and checks equality to every saved ambient Gram. Rational LDL decomposition verifies positive definiteness and determinant 1. It checks every normal Gram, saturation, every cluster distance, and exact rational squares of densities and comparison ratios.

Run:

```sh
python research/novel_routes_20261113/higher_parent/verify_gamma72.py
```

The replay requires only NumPy and SymPy. Its output is `gamma72_verified.json`.

## The comparison ceiling and its scope

If a rank-`k` normal lattice has minimum at least 8, the exact Hermite constant in dimensions `k<=8` gives

\[
\det K\ge(8/\gamma_k)^k.
\]

Consequently, a bare-section lower bound obtained by retaining radius `sqrt(8)/2` is at most

\[
2^{(72-k)/2}(\gamma_k/8)^{k/2}
=2^{36-2k}\gamma_k^{k/2}.
\]

For dimensions 71 down to 64 these ceilings are respectively

\[
2^{34},\quad 2^{33}/\sqrt3,\quad 2^{30}\sqrt2,\quad
2^{29},\quad 2^{27}\sqrt2,\quad 2^{27}/\sqrt3,\quad
2^{25},\quad 2^{24}.
\]

The new finite clusters beat these particular ceilings. This statement does **not** rule out a section whose actual minimum exceeds 8, a packing from a different parent, or another periodic construction.

## Discovery and reproducibility

- `simplex_gram_search.py` searches a finite family of negative-off-diagonal Gram matrices. Diagonals are 8 for ranks 2–5; off-diagonals are in `{-4,...,0}` for ranks 2–4 and `{-3,...,0}` for rank 5. The recorded 1,048,576 rank-5 matrices are a finite search, not a classification of all possible normal lattices. The short-vector screen in this abstract search is only a discovery filter.
- `embed_simplex.py`, `embed_grown_simplex.py`, and `embed_grown_simplex6.py` realize selected abstract Grams in the actual parent. Realization in the minimum-8 parent supplies the missing minimum guarantee for the normal lattice.
- `grow_pool.py` starts from basis vectors and short sums/differences and produces a capped pool of 50,000 antipodal minimum-8 rays. Every returned vector is checked exactly. This is a sample, not the full shell.
- `embed_roots.py` and `embed_grown_root.py` find primitive `4D5`, `4E6`, `4E7`, and `4E8` normal lattices.
- `root_clusters.py` discovers finite dual clusters using floating enumeration followed by exact integer norm/diameter calculations. Its greedy clique search is not a global maximum-clique proof. The final verifier needs none of the search completeness assumptions.

`normal_mutation.py` additionally attempted 60,000 integer normal-Gram mutations for each of the fixed nine-point 66D and 65D clusters (78,831 total passed its finite minimum screen). Neither improved on its starting value. This is a heuristic result, not a nonexistence claim.

Large reproducible discovery pools can be regenerated and are excluded from the checkpoint. The compact final witnesses are retained.

## Sources and priority investigation

1. Gabriele Nebe, *An even unimodular 72-dimensional lattice of minimum 8*, https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/papers/Ext72.pdf . This proves the inherited parent minimum. The primary paper is linked; its full PDF and extracted text are excluded from this checkpoint.
2. Nebe–Sloane catalogue, *The Lattice Gamma72*, https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/neb72.html . This supplies the exact Gram; archived as `neb72.html`.
3. Ruitao Chen et al., *New Sphere Packings from the Antipode Construction*, arXiv:2505.02394, https://arxiv.org/abs/2505.02394 . The method of suboptimal normal sections is established there for dimensions 19, 20, 21, 23, 44, 45, 47. Their stated results do not cover 65–71.
4. Henry Cohn, sphere-packing density table, https://cohn.mit.edu/sphere-packing/ . Its stated dimension coverage is 1–48, 56, 64, 72; it does not provide 65–71 record benchmarks.
5. Nebe–Sloane densest-packing table, https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/density.html . Above 48 its entries are sparse: 54, 56, 64, 72, 80, 128. It establishes the 64D control benchmark but cannot establish priority for the intervening dimensions.

Targeted searches for 65–71 density results and Gamma72 antipode applications found no direct prior comparison value. Search absence is not a proof of novelty. In particular, the benchmark investigation has not classified possible laminations of Ne64, higher-minimum sections, or sections of the Mordell–Weil 80-dimensional lattice.

## Earlier existential 65D benchmark from Ne64

Before the explicit deep-coset search, a baseline check found a stronger **existential lattice construction** in dimension 65. This invalidates any interpretation of the Gamma72 65D row as a record candidate.

Let `L=Ne64`, with minimum 12, determinant `3^32`, and `sqrt(3) L*` isometric to `L`. The standard theorem on theta series of strongly 3-modular lattices writes its theta series in the basis

\[
\theta_{A_2}^{32-6j}\bigl(\eta(\tau)\eta(3\tau)\bigr)^{6j},\qquad 0\le j\le5.
\]

Vanishing through norm 10 gives coefficients
`[1,-192,10944,-185088,536256,-483840]`. Exact power-series arithmetic gives counts

- norm 12 in `L`, equivalently norm 4 in `L*`: `138458880`;
- norm 18 in `L`, equivalently norm 6 in `L*`: `37527682805760`.

The discriminant group is a nondegenerate 32-dimensional quadratic space over `F3`. Regardless of its Witt type, it contains at least

\[
3^{31}-2\cdot3^{15}-1=617673367586132
\]

nonzero isotropic classes. Such classes contain only vectors of even norm. Since `L*` has minimum 4, the classes admitting a vector of norm below 8 are occupied by vectors of norm 4 or 6. Their number is at most `37527821264640`, strictly smaller than the available isotropic classes. Thus some nonzero class `h+L` has minimum at least 8. Its negative has the same minimum, and `3h` belongs to `L`.

The 65-dimensional lattice generated by `(L,0)` and `(h,2)` has minimum 12: nonzero layers modulo 3 have horizontal squared length at least 8 and vertical squared length at least 4; layers divisible by 3 have vertical squared length at least 36. Its covolume is `2*3^16`. Therefore

\[
\delta_{65}\ge\frac{3^{33/2}}2
=37,279,553.9356211\ldots,
\]

which is larger than the Gamma72 65D value. This proof depends on the cited published Ne64 properties and the standard 3-modular theta-series structure theorem, in addition to the exact arithmetic in `ne64_theta.py` / `ne64_theta.json`. It is an existence proof and does not select an explicit class representative `h`.

The primary catalogue source for Ne64 is archived as `CQ64b.html`, from https://www.math.rwth-aachen.de/homes/Gabriele.Nebe/LATTICES/CQ64b.html . The theta-series structure theorem is due to H.-G. Quebbemann, *Modular lattices in Euclidean spaces*, Journal of Number Theory **54** (1995), 190–202. The priority of this natural 65D extension has not been established.

The Ne64 formal theta coefficients also agree with the independently maintained Sloane sequence https://oeis.org/A034642 (checked in this run).
