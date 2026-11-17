# Independently audited deep-coset parity extensions in83–85D

Run `python run_deep.py` (numpy and python-flint) for exact source maps, full Gram determinant, parity geometry and outward-rounded density. `python run_deep.py --full` additionally compiles C++17 and replays inherited source and fresh depth exclusions.

The vertical lattices are A3, D4, D5 in dimensions83,84,85, respectively, all with Gram determinant4. Every vertical root parity has a horizontal class isometric to the previously independently certified source101 of depth at least6. Each remaining nonzero parity has integral norms and a fresh exhaustive exclusion of norms at most3, hence depth at least4. Zero parity has vertical minimum8. Exact parity enumeration and the inverse-Gram coordinate bound certify that no smaller vertical vector is missed. Together these imply squared minimum at least8. The pinned80D parent's minimum8 remains inherited from Bachoc–Nebe's theorem with its constructive hypotheses checked.

Full Grams and exact squared center densities are in `deep_extension/verified83.json`, `verified84.json`, `verified85.json`. Independent geometry, source mapping, local shell coverage and density checks are in `independent_audit/`. Global decoder searches were reviewed and freshly rerun, not independently reimplemented. Original historical PENDING_AUDIT labels are preserved; the independent report records completed checks.

| Dimension | Ordinary density at least | Derived BCS88 section | Gain over that comparison |
| --- | ---: | --- | ---: |
|83|3.08087078338732e-18|√2D5|>465.685424%|
|84|1.18808157306620e-18|√2D4|300% exactly|
|85|4.55474871517601e-19|√2A3|>182.842712%|

[83D source audit](../bcs83_comparator_20261117/); [84/85 source derivation](../priority_resolution_20261115/CLASSICAL_SECTIONS.md). These are derived classical constructions, not printed dimension-specific historical records or an exhaustive worldwide-priority comparison.

The included old run_all.py checks the inherited82D dependency package. The83D primary proof is in deep_extension/PROOF83.md; D4/D5 use the same horizontal-coset/vertical-parity argument with their explicitly stored Gram forms. Portable candidates retain the audited generators/known-class witnesses but remove unused auxiliary orbit entries; the original frozen candidate snapshots and hashes are retained in independent_audit/.
