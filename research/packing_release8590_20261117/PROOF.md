# Correlated Mordell-Weil dual-normal sections with optimized simplex translate clusters

This package contains six periodic sphere-packing constructions in dimensions 85–90, not upper bounds or optimality claims. It retains only the best independently audited input in each dimension. Comparisons are to the named, pinned repository certificates, not an exhaustive current-record survey.

## Inherited parent and exact reconstruction

The parent is the rank-92 Mordell-Weil lattice associated with `y²=x³+1+t^48` in characteristic 47. Its squared minimum is at least 16 and its Gram determinant is `47^14`. The minimum theorem is inherited from Schütt–Shioda, *Mordell-Weil Lattices*, Theorems 6.23 and 13.35 and §13.4.2; it is not proved by a heuristic vector search. Source: https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf .

The `parent()` routine retained in `packing_resume_20261117/verify.py` checks all 180 explicit polynomial sections against the curve, replays 16,110 pairwise heights by polynomial gcd, verifies the resulting Gram matrix, and checks that the coordinates generate the full integer module. Both verification routes rerun this routine. They share this parent reconstruction; the separate independence lies in the section saturation, coarea and density checks, not in an independently formalized parent theorem.

## Antipode construction

Use row coordinates for the parent Gram matrix `G`. Let `B` be the full-row-rank integer normal matrix, of size `k×92`, and let `S` be the finite set of integer labels in the frozen file. Smith/Hermite checks prove that the map `z ↦ z B^T` is onto `Z^k`.

Set `K=B G^(-1) B^T`, and let `Λ` be the integer kernel of `B`. Choose an integer lift of each label `s`, then orthogonally project it into the kernel subspace. This gives a translate `p_s+Λ`. All lifts and rational shifts are explicitly reconstructed by the primary verifier.

Define `β=max_(s≠t) (s−t) K^(-1) (s−t)^T`. Distinct points from two translates differ by a nonzero parent vector, whose squared length is at least 16. Its normal component has squared length at most β. Therefore the projected squared separation is at least `16−β`. Within a translate it is at least 16. All translates are distinct modulo Λ: otherwise a nonzero parent vector would lie entirely in the normal subspace with squared norm at most β<16, a contradiction.

The period Gram determinant is `det(G) det(K)`. Hence, for `n=92−k` and `M=|S|`, the exact center-density square is

`δ² = M² ((16−β)/4)^n / (47^14 det(K))`.

Occupied-volume density is `V_n δ`. The verifiers use rational Machin bounds for π and explicitly check the displayed lower floors. No numerical optimizer output is treated as a certificate.

## Verification routes and scope

The primary route constructs a unimodular Hermite kernel, integer label lifts, period Gram and rational projected shifts. The independent route checks Smith saturation, all exact pair losses and the coarea determinant formula. Both verify the exact density-square identity and conservative decimal lower rounding.

The input files and hashes are unchanged from the independent audits. The portable primary wrapper now selects the baseline explicitly by dimension. The independent 182-normal script is restricted to the four current 86, 87, 88 and 90 inputs; its superseded 85 input is intentionally omitted. The unchanged 192-normal 85 and 200-normal 89 inputs have separate audits. None of the unproved or rejected 93-dimensional experiments is included.

## Pinned comparisons

`RESULTS.json` records exact input hashes, safe occupied-density floors, and gains. Baselines are pinned to repository commit `20b638f5d5b0bf00a2a4ce8d70c4e062517a502e`:

- 85: `published_adaptive_results.json`, adaptive exceptional-coset lattice construction.
- 86–89: `published_simplex_results.json`, preceding MW92 simplex antipodes.
- 90: `published_mw92_results.json`, preceding three-translate MW92 antipode.

These baseline files are retained exactly for reproducible comparison. The public source attribution and broader comparison survey remain those of the corresponding repository proof packages.
