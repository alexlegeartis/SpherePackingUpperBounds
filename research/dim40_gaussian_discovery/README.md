# 40D discovery and deterministic export companion

The separately frozen proof is `../../verifications/dim40/verify_primary.py`. This companion contains the original numerical candidates, exact source data, settings and a deterministic reconstruction of the final rational witness. It does not claim that a fresh optimizer run is bitwise reproducible. The theorem does not depend on optimizer success or sampled feasibility.

## Verified replay
Run `python replay_export.py` using Python 3.11+ and the dependencies in requirements.txt, without `-O`. This re-exports the angular rational data, reconstructs the physical and Fourier radial polynomials, proves both complete radial rays and pair coupling, and compares all resulting JSON fields with the frozen certificate. It writes only `work/`. The complete cyclic angular domain must still be replayed by the proof package.

## Discovery settings
The source is the exact 44D radial degree-111 witness in inputs/source44, restricted to physical dimension 40 by the exact Fourier marginal identity H_40[j] = sum_{k>=j} H_44[k] binom(k,j) (2)_(k-j). The radial LP combines a nonnegative multiple of that source with 32 signed Laguerre coefficients (degrees 0 through 31), normalized Fourier value 1, nonnegative origin cost, squared minimum distance 21/4, physical cutoff 7, and kappa = 17/7.

The dimension-specific cyclic angular search uses degree at most 19, alpha = 40, nonnegative normalized Gegenbauer coefficients, c0 = 1, and eight cut rounds. Exact export rationalizes coefficients to denominator at most 10^12, subtracts 1/1000 from c0, then renormalizes to c0 = 1. The final polynomial has degree 18.

The radial LP uses normalized row margin 10^-8, HiGHS primal and dual feasibility tolerances 10^-9, initial fixed physical/Fourier grids and exact-root-guided added points. It rationalizes the final solution to denominators at most 10^14. The frozen 40D origin M receives a factor 5001/5000 (0.02%) safety cushion, distinct from the 0.2% default used in later dimensional experiments. Eight repair rounds preceded the recorded passing export. Intermediate rejection JSON is historical evidence of roots in earlier candidates, not a claim against the final witness.

## Optional fresh searches
`python search_angular.py --dimension 40 --degree 19 --alpha 40 --rounds 8` writes an exploratory angular candidate in work/.
`python search_radial.py --dimension 40 --degree 31 --r2 5.25` searches using the frozen angular polynomial. To seed the historical radial separation points, copy inputs/radial_scaled40_degree31_r5p25_extra.json to work/ first. Search outputs are always numerical candidates, never certified bounds. Floating-point solver/library versions can change their coefficients or success. No new optimizer replay was performed in preparing this companion.

The exact replay uses the candidate stored in inputs, not any optional new search output. See the separate proof package for the strict density claim 0.000007516663425 and its pinned comparator.
