# Sections, normal subspaces, and independent audits

## Status

**No new sphere-packing density record was found in this directory.** The searches produced exact finite obstructions, independently checked parameter calculations, reusable vector pools, and numerical evidence against several local deformation routes. A failed bounded search is not a nonexistence theorem.

All densities below are **center densities**: a lattice with squared minimum μ and Gram determinant d has center density μ^(n/2)/(2^n sqrt(d)). Multiplication by the unit-ball volume gives packing density. The previously certified 36D lattice in `constructions/dim36` remains the starting result; this directory does not claim that result as a new discovery.

The P48 ambient minimum6 and unimodularity are inputs from the published constructions and the external certificate repository. They were not independently re-proved by a complete 48D enumeration here. The external repository was read at commit `a0bc7f51444e42c3cbae0d2a910039cd50013500`:

- <https://github.com/wcgbg/sphere-packing>
- <https://arxiv.org/html/2607.20359v4>

The Kappa11 source is the Nebe–Sloane lattice catalogue:
<https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/KAPPA11.html>.

## Exact finite results

### A 23-contact obstruction to two binary parity checks

`mod2_obstruction_certificate.json` contains 23 explicit integer-coordinate vectors in the certified 36D lattice. Each has squared norm8 under `Q_original.txt /20`, and

\[
\sum_{v\in S}vv^T=0\pmod 2,\qquad |S|=23.
\]

Suppose two binary linear forms a,b were nonzero jointly on every listed contact. Over F2, the polynomial

\[
q(v)=a(v)+b(v)+a(v)b(v)
\]

would equal1 on all 23 vectors. The displayed moment identity makes the sum of every degree-one or degree-two monomial zero; diagonal entries also give the degree-one identity. Thus the sum of q would be both0 and1. This is a short, exact contradiction.

The conclusion concerns maps to F2² and sublattices defined by at most two binary parity checks. It does **not** rule out integer hyperplane sections, cyclic quotients of order4, or higher-rank parity constructions. No completeness assertion about the contact shell is needed for this proof.

### Independent audit of the fixed length17 mixed-code obstruction

`audit_n17.py` independently reconstructs the code in `../mixed_codes/n17k9_G.npy` using polynomial multiplication modulo x²+x+1. It does not use the main C++ enumerator or its output to construct the code. Its 262,144 codewords have minimum weight7, with 1,224 words of weight7 and 408 distinct weight7 supports. The support set agrees exactly with the archived set and is a 2-(17,7,63) design; each coordinate occurs168 times.

For nonzero parity columns in F3², exhaustive signed-sum arithmetic shows that a seven-coordinate support avoids a full-support signed dependency exactly when its four projective color counts are a permutation of (6,1,0,0). Every support would then contain15 monochromatic coordinate pairs. Double-counting those pairs gives

\[
63M=408\cdot15,\qquad M=680/7,
\]

contradicting integral M. This excludes rank-at-most-two parity matrices with nonzero columns for this particular outer code. The no-zero-column condition matters: in the intended mixed hexagonal construction a zero parity column permits a ternary weight-one word with squared cost4, below the target minimum8. The statement does not cover arbitrary length17 outer codes or rank3 parity matrices.

`audit_n17.json` records the independent weight enumerator, source hash, design parameters, all12 ordered admissible color-count patterns, and the contradiction.

### Nearest integral roundings of scaled Kappa11

The ideal scaled form (3/2)Kappa11 has minimum6 and determinant84,075.75. That determinant would be small enough for a record-beating 37D P48 normal lattice **if** an appropriate primitive embedding existed. Its half-integral off-diagonal entries prevent the proposed direct integral embedding.

There are14 odd off-diagonal entries in the catalogue Gram. All 2^14 =16,384 nearest integral roundings were checked against216 explicit exact norm-four Kappa11 vectors. Every rounding makes at least one of those vectors shorter than6. This is an exact finite obstruction to those nearest roundings, not a classification of all integral rank11 forms.

Wider integer MILPs gave positive-definite candidates of determinants145,200,134,784,139,920. A floating 11D discovery enumeration found no vector of norm below6 in those candidates. None meets the determinant target102,382.706675..., and no P48 embedding was found or asserted. A log-determinant iteration stalled at134,784 and was stopped; its repeated identical iterates are not independent search starts. The script now stops after16 unchanged iterates.

### Exact replay

From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/density_breakthroughs_20261113/sections/verify_finite_results.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/density_breakthroughs_20261113/sections/audit_n17.py
```

These commands replay exact witness norms, the binary moment obstruction, nearest-rounding obstruction, saved normal determinants and saturation indices, displayed antipode cluster distances, and density arithmetic. They deliberately do not turn floating enumeration or numerical linear algebra into a proof. The first writes `finite_verification.json` and the second writes `audit_n17.json`.

## Normal-subspace searches in P48

For a primitive rank-k normal lattice K in a unimodular P48 lattice, its orthogonal section has dimension n=48-k and determinant detK. The guaranteed bare section density is

\[
\delta_n=(3/2)^{n/2}/\sqrt{\det K}.
\]

The best saved normal spaces across these searches are:

| Section dimension | Normal rank | Best exact detK | Center density from min6 | Benchmark | Record-beating determinant threshold |
|---:|---:|---:|---:|---:|---:|
| 37 | 11 | 121,500 | 5.192782752 | 5.656854249 | less than102,382.706675 |
| 36 | 12 | 254,016 | 2.932325159 | 5.281629001... | not competitive |
| 35 | 13 | 391,500 | 1.928551098 | 2.828427125 | less than182,013.700756 |
| 34 | 14 | 648,675 | 1.223313656 | 2.222034570 | less than196,608 |

Every saved best basis was checked by exact Gram arithmetic and Smith normal form. The displayed saturation indices are1. The scripts also divide by the squared index when an unsaturated basis occurs. The search objective uses floating Schur complements, while reported improved determinants are recomputed exactly.

The search coverage was:

| Run | Explicit pool | Moves and scope | Completed budget |
|---|---|---|---|
| `p48_sections.py` | 8,192 P48p norm6 representatives, generated from published normal bases and exact root exchanges | Random starts; one-, two-, and three-row deletion/refill exchanges; ranks11–14 | 1,481 trials,300 seconds |
| `p48_large_sections.py` | 172,204 distinct P48p representatives, merging the small pool and a capped shell | Starts from all or five published K10 rows, or one random pool vector; one-, two-, and three-row exchanges | 263 trials,300 seconds |
| `p48n_extended_sections.py` | 600,454 P48n vectors:300,457 of norm6;149,997 of norm8;120,000 of norm10;30,000 of norm12 | Starts from all or partial E8 bases or a random vector; one- through four-row exchanges; ranks11,13,14 | 36 trials,180 seconds |

The P48p shell discovery stopped at its300,000,000-node cap, producing167,096 antipodal representatives in27.7 seconds. All returned leaves were individually rechecked to have exact norm6. It is an incomplete pool. The P48n pool was supplied by the neighboring `gluing` workstream and has individually exact norms and ambient integer coordinates. Neither pool is asserted to contain the full P48 minimum shell.

The published K10-containing rank11 extension floor is121,500. The later searches deliberately allow deletion of the K10/E8 rows and arbitrary new normal subspaces; their failures are not extensions of that published nonexistence result.

### Rank11 dual clusters

`K11_extension.npz` stores an explicit rank11 normal lattice of determinant121,500. A floating discovery run of its dual through squared norm1 returned545 antipodal representatives in7,050 nodes without hitting the cap. The output vectors and all graph distances are exact; floating pruning was not converted to an outward interval completeness certificate.

`dual_clusters.py K11` solves the finite origin-containing clique problems for11 attained diameter thresholds from9/20 to1. Every finite graph search ended uncapped, with at most119,711 recursive calls. The largest cluster at diameter1 has12 points. The best non-singleton density among the tested threshold optima is2.957648570, from four points at squared diameter3/5; it is worse than the bare section density5.192782752. `verify_finite_results.py` checks the displayed clusters directly, without claiming optimality over the infinite dual lattice.

For perspective, a two-point antipode cluster with primitive dual difference y and squared diameter β induces the lower-rank normal lattice K0=K intersect y-perp with detK0=detK·β. Its nD density divided by the bare (n+1)D section density is

\[
4\sqrt{\beta/6}(1-\beta/6)^{n/2}
\le \frac4{\sqrt{n+1}}\left(\frac n{n+1}\right)^{n/2}.
\]

A nonprimitive dual difference only decreases this ratio. Thus a two-point gain in these dimensions would require a substantial gain in an induced section one dimension higher. This does not exclude clusters with higher affine dimension.

## Numerical screens around the 36D seed

Six principal 35D sections, deleting original coordinates20,29,30,33,34,35 (zero based), have minimum8 and center density2.470255771. Their dual normals have squared norm16/7. Explicit contacts from the shared floating norm8 shell give full numerical rank630 in the symmetric35D form space. Positive numerical eutaxy weights were found with minimum weight at least1.2680e-4 and residual at most7.3e-15. These are strong local numerical screens, not exact eutaxy certificates.

The shared source shell has94,050 antipodal representatives. Individual contact norms are exact, but this directory does not claim an independent proof of shell completeness. `contact_screen.py` records an earlier incomplete-contact LP screen that timed out; the later `eutaxy_screen.py` results supersede it.

For a hypothetical minimum10 section of the new36D lattice, a primitive 35D dual normal could beat the35D benchmark if its squared norm were less than4297.998464583. A primitive codimension-two normal lattice could beat the34D benchmark if its determinant were less than2785.573433035. These are conditional density thresholds, not constructed packings.

`normal_section_search.py` searches integer normal pairs in a saved dual basis, minimizing the number of sampled norm8 contacts in their common kernel. A90-second soft-penalty search completed7 starts; a60-second hard-cap search completed34 starts. The best hard-cap pair still retained1,065 norm8 lines at normal determinant2778.883786848. No pair passed even the necessary norm8 filter. Avoiding all norm8 contacts would still leave the norm9 shell to exclude before claiming minimum10. The recorded historical cap is approximately1e-6 more restrictive than the exact threshold; the script now uses the exact-density decimal threshold.

## Full 43D periodic screen and normalization audit

`periodic43_fista.py` independently reconstructs the six-coset P48n packing used in `../contact_walk/periodic43.py`. Exact checks confirm that the normal basis is saturated and has determinant3072. The active normal squared differences are7/24, so the contact squared distance is137/24 and center density is

\[
6(137/96)^{43/2}/\sqrt{3072}=226.517274889\ldots.
\]

The variable space contains946 symmetric metric entries and215 relative coset displacements, totaling1,161. Global translation is fixed by holding coset0 stationary. The contact derivative was independently checked to be

\[
v^TDv+2v\mathbin{\cdot}(s_i-s_j).
\]

The objective trace(D) has the correct sign: at fixed minimum distance, the density's logarithmic derivative is -trace(D)/2. The Cholesky and null-space orientation, off-diagonal factor2, coset-pair signs, and ambient-to-section distance subtraction were audited.

The ambient pool supplies75,964 exact contact vectors over12 active coset pairs. A fixed-seed sample of12,000 contacts has numerical Jacobian rank1,161. Preconditioned FISTA found nonnegative weights after700 iterations, with11,985 entries positive and Euclidean residual

\[
\|A^T w-\mathrm{trace}\|_2=1.35\times10^{-11}.
\]

This is numerical evidence against first-order density descent in the full metric-plus-translation space. It is not an exact local-optimality proof and does not exclude nonlinear finite deformations. No claim of complete contact enumeration is needed to interpret a valid finite positive-weight certificate, but the displayed weights remain numerical.

The compact `periodic43_fista.npz` is412,318 bytes. It stores only the12,000 selected exact ambient coefficient rows (int16; maximum absolute coefficient17), pair labels, nonnegative weights, and G/T/E reconstruction data. It contains no dense contact Jacobian or ambient pool. A self-contained replay, requiring no ignored pool or external clone, is:

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python research/density_breakthroughs_20261113/sections/check_periodic43.py
```

This reconstructs the Jacobian, verifies individual ambient norms and oriented normal projections exactly, and rechecks the numerical weight residual and rank. Its output is `periodic43_witness_checked.json`.

## Reproduction of exploratory runs

Python dependencies are NumPy, SciPy, SymPy; `p48_enumerate.py` additionally needs fpylll and cysignals. The scripts use fixed random seeds. Wall-clock-limited runs can complete different numbers of trials on different machines; replay of saved finite witnesses is deterministic.

The following commands assume the external repository and neighboring gluing inputs are available. On the original shared machine these jobs were pinned to CPUs2–3, with at most two BLAS threads.

```bash
g++ -O3 -std=c++17 research/density_breakthroughs_20261113/contact_walk/shell.cpp -o ../shell_discovery
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/density_breakthroughs_20261113/sections/p48_pool.py P48p 8192
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/density_breakthroughs_20261113/sections/p48_enumerate.py P48p
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/density_breakthroughs_20261113/sections/p48_sections.py P48p 300
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/density_breakthroughs_20261113/sections/p48_large_sections.py 300
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/density_breakthroughs_20261113/sections/p48n_extended_sections.py 180
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/density_breakthroughs_20261113/sections/prepare_section_inputs.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python research/density_breakthroughs_20261113/sections/normal_section_search.py 60 81723 hard
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python research/density_breakthroughs_20261113/sections/periodic43_fista.py 12000
```

`prepare_section_inputs.py` requires the shared `contact_walk/contacts_reduced.txt` shell and the compact saved `dual_lll.npz`. The root `contact_walk/prepare_shell.py` and floating shell executable regenerate that discovery shell. Do not interpret that regeneration as a new minimum certificate.

## Compact deliverables and regenerable files

Retain this README, `audit_root_claims.md`, Python source files, JSON reports, `mod2_obstruction_certificate.json`, `audit_n17.json`, the four best normal-basis files used by `verify_finite_results.py`, `K11_extension.npz`, `K11_clusters.json`, `kappa11.json`, `kappa11_vectors.txt`, and `dual_lll.npz`. Together these preserve the exact finite claims and the search conclusions. The compact `periodic43_fista.npz` and its self-contained checker preserve the numerical contact/weight witness and should be retained.

Large discovery pools and duplicate intermediate arrays can be ignored in a compact checkpoint: `P48p_enumeration/vectors.txt`, `P48p_enumeration/pool.npz`, `P48p_pool.npz`, `full_contacts.npz`, `normal_contact_coordinates.npz`, and any regenerated `eutaxy_*.npz`. Their generation code and finite statistics are retained. The neighboring external repository and P48n gluing pools are also regenerable inputs rather than new certificates. No registry, existing construction certificate, or benchmark table was changed by this work.
