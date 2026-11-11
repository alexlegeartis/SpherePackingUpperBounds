# 38D discovery and deterministic export companion

The separately frozen proof is `../../verifications/dim38/verify_primary.py`. This companion contains the numerical candidates, exact source data, settings and a deterministic reconstruction of the final rational witness. No fresh optimizer replay was performed in preparing it. Optional fresh searches may change coefficients or fail because of numerical/library differences; they are never proofs.

Run `python replay_export.py` using Python 3.11+ and requirements.txt, without `-O`. It reconstructs angular rational data and radial physical/Fourier polynomials from frozen inputs, verifies both radial rays and pair coupling, and asserts complete JSON equality with the frozen certificate. It writes work/. Complete angular cyclic coverage remains in the proof package.

## Discovery settings
The exact 44D source polynomial (degree 111) is restricted to physical dimension 38, with Fourier marginal coefficients H_n[j] = sum_(k>=j) H_44[k] binom(k,j) ((44-n)/2)_(k-j). The radial LP adds a signed normalized Laguerre correction of degree at most 31 and optimizes a nonnegative source multiplier and nonnegative origin cost. The Fourier value at zero is fixed to 1. Squared minimum distance s=19/4, physical cutoff 7, kappa=51/19, and kappa*pi*s=51*pi/4>40.

The angular LP uses degree at most 19, alpha 40, nonnegative normalized dimension-38 Gegenbauer coefficients and c0=1. Export rationalizes to denominators at most 10^12, subtracts 1/1000 from c0, and renormalizes. The radial final search uses normalized row margins 10^-6, HiGHS method highs and primal/dual tolerances 10^-9, fixed initial grids plus exact-root-guided separation points. Earlier searches sometimes used the weaker 10^-8 margin; historical rejects are not results. Radial coefficients are rationalized to denominators at most 10^14, and M receives the exact factor 501/500 (0.2%).

## Optional fresh searches
`python search_angular.py --dimension 38 --degree 19 --alpha 40 --rounds 8` writes a numerical angular candidate in work/.
`python search_radial.py --dimension 38 --degree 31 --r2 4.75` uses the frozen angular witness. To seed historical root separation points, copy inputs/radial_scaled38_degree31_r4p75_extra.json into work/ first. No fresh optimizer replay is implied.

The replay uses frozen inputs, not optional new search output. It is not a second independent mathematical audit.
