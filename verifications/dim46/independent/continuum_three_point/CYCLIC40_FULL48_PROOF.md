# Optimized cyclic 48D certificate: objective below17.835

Producer verification and independent audit both pass. This supersedes the intermediate18.779 candidate, and is stronger than the independently audited18.912 certificate.

The complete construction is CYCLIC_FULL48_PROOF.md with these explicit changes:

- Angular polynomial: cyclic_angular40_exact.json, degree19, c=5001/10000, B approximately907537125.6591902.
- It has nonnegative normalized48D Gegenbauer coefficients, constant coefficient1, is less than-10^-5 on[-1,1/2], negative on[-1,c], and has exactly one root on(c,1) before becoming positive. rationalize_cyclic_angular40.py verifies these facts by rational polynomial construction and Sturm counts, from the fixed discovery snapshot cyclic_angular40_discovery_snapshot.json.
- The single-sign-change scale lemma in CYCLIC36_PROOF.md applies unchanged. At alpha40 the actual triangle expression is negative over the compact domain, verified by certify_cyclic_optimized_actual40.py with argument40:11,932 rectangles, depth19,128-bit precision. For t<=0 and u>=5 the same analytic tail applies, with c=.5001 giving first gap>4.001 in the latter case. The exact exponential estimate40^30/30!>2B*10^5 is sufficient. Thus the global cyclic sign holds for every alpha>=40.
- Radial data: adaptive_cyclic40_radial111.json, degree111, normalized by its Fourier constant coefficient.
- kappa=17/8, h(x)=exp(-(17/8)pi(|x|²-6)). After scaling the minimum distance to1, alpha=51pi/4>40.
- M is neighbor_cost/(B H(0)), multiplied by1000001/1000000, as recorded exactly in cyclic40_radial111_exact.json.
- The same PSD feature construction gives f3(0,0)=MB and the full coupled condition f2<=M(2h-h²).
- certify_cyclic40_radial111.py verifies both global polynomial signs and the coupled interval6<=|x|²<=7 by33 rational/Taylor/Arb subdivisions at384 bits. Beyond7 the radial polynomial is negative and0<h<=1, so coupling follows automatically.

The exact normalized objective is17.834291498319086... <17.835. Multiplication by pi^24(3/2)^24/24! gives the packing density upper bound. The existing continuum theorem and sign proof establish a density bound only; no lattice construction, optimality, kissing-number conclusion, or novelty/priority claim is implied.
