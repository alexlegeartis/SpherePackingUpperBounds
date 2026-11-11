# Comparison sources and conventions

Reviewed 10 November 2026. Except for the updated 4D and 5D rows below, the numerical baseline is [Henry Cohn's sphere-packing table](https://cohn.mit.edu/sphere-packing/), retrieved through the search index (reported crawl age: four months). Direct live retrieval failed with HTTP 502. Thus these are pinned, source-attributed comparisons, not a claim that no newer results exist.

## Lower bounds

- **4D:** the D4 lattice has density pi^2/16, as tabulated by [Cohn](https://cohn.mit.edu/sphere-packing/). This is the known construction, not an optimality claim for unrestricted packings.

- **21D:** Chen, Hu, Li, Wang and Wu, [New Sphere Packings from the Antipode Construction](https://arxiv.org/abs/2505.02394v1), 2025. Center density 2^(-41)*3^(-23/2)*43^(21/2), as listed in Cohn's general-packing table. Multiply by 4^11*11!*pi^10/22! for density L.

- **19D:** Chen, Hu, Li, Wang and Wu, [New Sphere Packings from the Antipode Construction](https://arxiv.org/abs/2505.02394v1), 2025. Center density 3^(-9) * 5^(-21/2) * 13^(19/2), as listed in Cohn's general-packing table. Multiply by the unit-ball volume 4^10*10!*pi^9/20! for the density L. This is a cited construction, not a new lower bound.

- **20D:** Chen, Hu, Li, Wang and Wu, [New Sphere Packings from the Antipode Construction](https://arxiv.org/abs/2505.02394v1), 2025. Non-lattice periodic packing, five spheres per underlying lattice cell. Exact density L = (pi^10/10!) * 3^20/(2^10 * 5^(21/2)).
- **48D:** Leech and Sloane, [Sphere packings and error-correcting codes](https://doi.org/10.4153/CJM-1971-081-3), Canadian Journal of Mathematics 23 (1971), 718–745. Extremal lattice packing density L = (3*pi/2)^24/24!. The attribution and center density are independently tabulated by Cohn (reference 21). Lattice constructions also supply lower bounds for unrestricted packings; they do not establish optimality among all packings.

## Previous upper bounds

**4D current comparison:** [arXiv:2512.04829v2, Table 1](https://arxiv.org/html/2512.04829v2#S3.T1), version dated 8 December 2025, reports 0.63610277 for unrestricted packing density. Checked directly on 10 November 2026. This is the paper’s printed value at eight decimal places; no additional digits, directed-rounding guarantee or independent certificate replay are asserted. Percentages therefore compare against that reported decimal, not a newly certified exact baseline. The repository bound <0.635766243 remains smaller. This paper reports dimensions 4–7 and 9–16 in Table 1, with a separate 8D case study; it provides no replacement comparison for our 19D, 20D, 21D or 48D rows.

**4D historical comparison:** Cohn–de Laat–Salmon, [Three-point bounds for sphere packing](https://arxiv.org/abs/2206.15373), supplies the earlier upper 0.6361073321551329 from Cohn's table. The new package also assumes four published near-neighbor average-degree bounds from [de Laat–de Oliveira Filho–Vallentin, Table 2](https://arxiv.org/pdf/1206.2608). See the 4D proof for exact thresholds and the limits of its independent computational audit.

For dimensions 19, 20, 21 and 48, Cohn's table attributes its upper bounds to Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini, [High-dimensional sphere packing and the modular bootstrap](https://arxiv.org/abs/2006.02560v2), JHEP 12 (2020), 066, [DOI](https://doi.org/10.1007/JHEP12(2020)066). The table's upward-rounded entries are 0.0098179551395438 (19D), 0.0071270536033763 (20D), 0.0051596603948176 (21D), and 0.0000004529067791 (48D). These concern unrestricted packing density, not a lattice-only bound. We use these published-table decimals, not the unknown exact optimum of the Cohn–Elkies optimization problem.

## Arithmetic and interpretation

Use the same lower bound L before and after the update, the referenced older upper U, and the repository's displayed certified strict upper N:

- Upper-bound decrease (%) = 100 * (U - N) / U.
- Additive-gap decrease (%) = 100 * ((U - L) - (N - L)) / (U - L) = 100 * (U - N) / (U - L).

The second percentage is the fraction of the earlier uncertainty interval removed. It is not a claim that N lies that percentage above the optimum. Both comparisons use densities, not center densities. N is rounded outward, so we do not take credit for unreported extra digits of the certificate.

`python comparison.py` checks all displayed values, exact-rational upper-percentage arithmetic, gap identities and agreement at 70 and 100 decimal digits. `python comparison.py --write` regenerates README/results tables and the registry. Lower formulas are numerical evaluations of cited exact expressions, not newly proved lower-bound constructions. The readable tables show densities to five significant digits with an approximation mark and percentages to two decimal places. These display values are not substitute rigorous bounds. The machine-readable registry retains its previous precision, including six decimal places for percentages; calculations use the unshortened inputs. `RESULTS.json` also records both additive gaps.

## Dimension5 addition

The lower density is pi²/(15sqrt(2)), from the standard D5 lattice, tabulated by [Cohn](https://cohn.mit.edu/sphere-packing/). The current comparison uses the paper-reported0.5122645 from [arXiv:2512.04829v2 Table 1](https://arxiv.org/html/2512.04829v2#S3.T1), not the older0.5126451306253027 three-point benchmark. We have not independently replayed the newer paper's own numerical certificates. The decimal is a pinned reported comparison, not an exact optimum claim.

The new5D certificate uses the rigorous neighbor-count results in [de Laat–de Oliveira Filho–Vallentin, Table2](https://arxiv.org/html/1206.2608#S6.T2) as external mathematical inputs. Their use is audited; those original spherical SDP certificates are not reconstructed here. The new witness is verified over complete continuous domains. Its audit is source review plus independent reruns and objective reconstruction, not a clean-room implementation or formal proof.

## Four-dimensional refinement

The previous p4 repository witness has outward upper 0.634526023. The earlier audited repository certificate 0.635766243 is preserved in the parent package. Relative to the updated paper-reported 0.63610277 baseline and D4 lower pi²/16, the conservative rounded comparison removes 8.189832% of the additive gap and lowers the upper by 0.247876%. The same four externally published neighbor-count inputs are required; audit consists of source review, isolated replay and independent objective reconstruction, not a clean-room implementation or formal verification.

The earlier 5D refinement uses the same external neighbor-count inputs and yields the outward-rounded bound0.510584101, closing3.5747932% of the gap to the D5 lower density relative to the pinned paper-reported0.5122645. Audit scope and benchmark caveats above remain unchanged.

The earlier radial-degree-six 4D witness strengthens the outward bound to **0.633430327**. Using this rounded value and the same updated paper-reported baseline, the additive-gap decrease is 13.881022%. Previous p3 and p4 certificates are preserved. The external count assumptions and audit scope are unchanged; the new far cosine cap 1/2 is proved by a radial-corner convexity argument in the refinement proof.

The earlier radial-degree-four 5D witness strengthens the outward bound to **0.510522525**, removing 3.705787% of the additive gap against the same paper-reported 0.5122645 baseline. Both previous 5D certificates are preserved. Neither the current 4D nor 5D result reaches a 20% gap decrease.

## Origin-coordinate refinements

The current registered bounds are 4D <0.633428734445, 5D <0.510501094478, 19D <0.009618343139, 20D <0.006979497272, 21D <0.005149549312. They result from exact Schur-complement slack removal in the existing PSD kernels, not new dimensions or a new theorem. Full source-domain replays and independent exact factorization and Arb objective reconstruction passed. Comparators and external mathematical inputs are unchanged; the generated tables recompute percentages from the strict outward decimals. See [proof and replay](verifications/origin_slack/).

## Dimension 18 addition

The certified upper is **<0.013408931**, independently replayed by two full-domain verification paths. Its comparator 0.0134853404450862 is the pinned Cohn indexed-table value attributed to [Afkhami-Jeddi et al.](https://arxiv.org/abs/2006.02560v2). The index was reported as four months old and live Cohn retrieval failed. Refreshed searches found no stronger unrestricted 18D bound, but no worldwide-priority claim is made. The newer arXiv:2512.04829v2 table covers other dimensions, not 18.

The lower density is pi^9/9! * (3/4)^9, from the center density (3/4)^9 of the 18D packing in Bierbrauer–Edel, [Dense sphere packings from new codes](https://doi.org/10.1023/A:1008733715001), whose abstract explicitly states the center density. The pinned Cohn index corroborates the construction (512 spheres per cell). This is a cited lower bound, not an optimality claim. The upper decrease is 0.566611% and the additive-gap decrease 1.044210%; comparison.py verifies the formulas.

## Dimension 17 cyclic addition

The strict upper **<0.018335552** is compared against the pinned Cohn indexed-table value 0.0184640903350649 attributed to [Afkhami-Jeddi et al.](https://arxiv.org/abs/2006.02560v2). Fresh searches on 11 November 2026 found no stronger finite unrestricted 17D bound, but the index was four months old and this is not worldwide-priority certification. The newer 2512.04829 table does not cover 17D.

The lower density is pi^(17/2)/(16 Gamma(19/2)), from the Lambda17 lattice's center density 1/16. This is tabulated by [Cohn](https://cohn.mit.edu/sphere-packing/) and [NumberDB](https://numberdb.org/Densities_of_the_densest_known_lattice_sphere_packings), with primary reference John Leech, [Some Sphere Packings in Higher Space](https://doi.org/10.4153/CJM-1964-065-1), Canadian Journal of Mathematics 16 (1964), 657–682. The construction is a lower bound, not an optimality claim. The generated table uses the same lower value on both sides.

The new proof uses a cyclic generalization: see [PROOF.md](verifications/dim17/PROOF.md) and [independent review](verifications/dim17/REVIEW.md), rather than claiming the literal individual-kernel-sign condition of Theorem 1.1. No neighbor-count input is used.

## Dimension 18 cyclic refinement

The current bound is **<0.013284728**, using the same lower construction and pinned published comparator as the original 18D entry. It lowers that comparator by 1.487634% and closes 2.741565% of the additive gap. The prior <0.013408931 certificate remains preserved. The cyclic generalization and all exact-domain checks are documented in [the refinement](verifications/dim18/refinement_cyclic/); no worldwide-priority claim is made.

## Dimensions 19 and 20 cyclic refinements

The current strict bounds are **19D <0.009557491** and **20D <0.006977702**, compared with the same cited lower constructions and pinned published upper values already documented above. Earlier repository bounds <0.009618343139 and <0.006979497272 are preserved in the origin-slack package. comparison.py recomputes every percentage and additive gap from the unshortened source inputs.

These improvements use the explicit cyclic generalization, not the pointwise-sign hypothesis of the printed Theorem 1.1. Full producer and independent exact-domain checks must pass. They add no neighbor-count assumption and make no worldwide-priority or optimality claim.

## Repaired 17D and cyclic origin saturation

Current strict bounds are 17D <0.018204096, 18D <0.013224473, 19D <0.009521000 and 20D <0.006936042. The 17D proof repairs the full Gaussian-polynomial physical tail using a second Gaussian scale, accounts for the changed Fourier normalization, and then saturates the isolated origin coordinate. The other three bounds use exact origin saturation of already-audited cyclic sources. None follows merely by substituting a solver objective.

The cited lower constructions and pinned published upper comparisons are unchanged. All percentages and gaps are regenerated by comparison.py; previous repository values remain in RESULTS.json and their complete proofs remain replayable. No optimality or worldwide-priority claim is made.

## Dimension 47 transfer and comparison

The strict certified upper is **0.000000616879891**. It comes from exact Fourier restriction of the 48D Schwartz witness, conversion to dimension-47 normalized Gegenbauer coefficients, and an exact rank-one origin Gram saturation. All inherited sign, coupling, cyclic-domain and analytic-tail proofs are replayed, then a separate independent rational implementation checks the transfer and rounding. This is an unrestricted packing upper bound, not a lattice-only result, optimality claim, external peer review or proof-assistant formalization.

The pinned published upper **0.0000006490757338** is from Cohn's indexed table citing [Afkhami-Jeddi et al.](https://arxiv.org/abs/2006.02560v2). On 11 November 2026 the indexed table was four months old and direct retrieval failed with HTTP 502. No exhaustive worldwide-priority claim is made.

The lower construction is Chen, Hu, Li, Wang and Wu, [New Sphere Packings from the Antipode Construction](https://arxiv.org/html/2505.02394v1), section “Dimension 47.” Its center density is 47^(47/2)/2^118. Thus the ordinary density is L = pi^23 * 47^23 * sqrt(47) * 24! / (2^70 * 48!), approximately 2.2669603229310984383e-8. This formula was checked against the primary construction and Cohn's indexed table; it is not a new lower-bound construction.

Using the outward upper, the upper reduction is **4.960260%**, and the additive-gap reduction is **5.139771%**. See [the proof](verifications/dim47/PROOF47.md).

## Dimensions 45 and 46

The certified strict bounds are 45D <0.000001311423575 and 46D <0.000000914882847. Pinned published upper comparators are 0.0000013300905665 and 0.0000009295151556, respectively, from the indexed Cohn table citing Afkhami-Jeddi et al., arXiv:2006.02560v2. The same table freshness and worldwide-priority caveats apply.

For 46D, Conway and Sloane, [The Antipode Construction for Sphere Packings](https://doi.org/10.1007/BF01232378), Inventiones Mathematicae 123 (1996), 309–313, supplies center density 13^23/3^(93/2). The ordinary lower density is (13*pi)^23/(23!*3^46*sqrt(3)), approximately 2.8609574721360648103e-8. The primary construction and indexed table agree. Upper reduction is 1.574187%, additive-gap reduction 1.624178%.

For 45D, [Chen, Hu, Li, Wang and Wu (2025)](https://arxiv.org/abs/2505.02394v1) supply center density 23^(45/2)/2^(183/2). The ordinary lower density is pi^22*23^22*sqrt(23)*23!/(2^45*sqrt(2)*46!), approximately 3.5582252942557527960e-8. The later [arXiv:2607.20359v4](https://arxiv.org/abs/2607.20359v4) reports a kissing-number improvement but explicitly ties this density, rather than improving it. Upper reduction is 1.403438%, additive-gap reduction 1.442014%.

46D retains the inherited angular cyclic proof and changes the radial Fourier pair. 45D also changes the angular polynomial and independently verifies its complete cyclic domain and analytic tails. Exact global radial sign proofs, coupling intervals, dimension-specific Gegenbauer positivity, normalization and origin Gram saturation are required in both dimensions. These are general-packing bounds, not lattice-only bounds.

## Dimension 6

The certified strict upper is 0.410189987. The baseline 0.4102905 is the printed general-packing value from [arXiv:2512.04829v2, Table 1](https://arxiv.org/html/2512.04829v2#S3.T1), not an independently reconstructed literature certificate. The E6 lattice gives ordinary lower density pi^3/(48*sqrt(3)), approximately 0.37294754558206493956. Thus the upper decrease is 0.024498%, and the additive-gap decrease is 0.269162%.

The count-enhanced cyclic proof uses published near-neighbor count inputs, cited in the package proof and [de Laat–de Oliveira Filho–Vallentin](https://arxiv.org/html/1206.2608#S6.T2). These external spherical SDP results are mathematical dependencies; their original certificates are not reconstructed here. The new certificate's exact PSD, Fourier identities, Gaussian guards, Schur factorization and complete pair/cyclic/far domains are replayed. The separate audit independently reconstructs algebra, guards and full pair inequalities and reviews the domain proof; it is not a proof-assistant formalization.

## Improved six-dimensional certificate

The current 6D strict upper is **0.409135399**, replacing the earlier registered 0.410189987 while preserving its proof. The benchmark remains the paper-reported 0.4102905 and the E6 lower remains pi^3/(48*sqrt(3)). Using the outward rounded current bound, the upper decreases **0.281532%**, and the additive gap decreases **3.093223%**.

The new frozen witness needs no Gaussian sign repair: eta = zeta = 0, and exact Sturm tests prove both full polynomial sign conditions. Exact origin Schur saturation, the cyclic continuum argument and the same external published near-neighbor count bounds remain required. The final witness arose from a 99/100 continuation mixture, but its proof checks the resulting exact coefficients directly; numerical feasibility of its parents is not assumed. The independent full replay checks 1,672 pair leaves, 4,359 far boxes and 87 compact cyclic boxes. The cited external count SDP certificates are not reconstructed. See [the new proof](verifications/dim6/refinement_mixture99/PROOF.md).

## Dimension 44 and small 6D follow-on

The strict 44D upper is **0.000001842868599**, compared with the pinned indexed-table upper 0.0000019017703144. The exact known lower from [Chen, Hu, Li, Wang and Wu (2025), Dimension 44](https://arxiv.org/html/2505.02394v1) has center density 157^22 * 2^-22 * 5^(-43/2) * 11^-23, hence ordinary density (157*pi)^22 / (22! * 2^22 * 5^21 * sqrt(5) * 11^23), approximately 3.9244265335277118067e-8. This was checked against the primary construction and the indexed Cohn row. The outward upper decreases the pinned upper by **3.097204%** and the additive gap by **3.162464%**. The upper comparison retains the table freshness and no-worldwide-priority caveats.

The current 6D follow-on upper is **0.409128533**, a small further gain over 0.409135399. Relative to the unchanged paper-reported 0.4102905 and E6 lower, the upper decrease is **0.283206%**, and the additive-gap decrease is **3.111610%**. Both Gaussian guard weights remain zero, and the external near-neighbor count dependencies remain unchanged. Earlier witnesses are preserved and replayed.

## Dimension 7

The strict upper is **0.320584772**, versus the printed value 0.32111475 in [arXiv:2512.04829v2, Table 1](https://arxiv.org/html/2512.04829v2#S3.T1). The standard E7 lattice lower is pi^3/105, approximately 0.29529787314571257310. The upper decreases **0.165043%**, and the additive gap decreases **2.052835%**. This is a comparison against a paper-reported decimal, not an independently reconstructed literature witness or exhaustive priority claim.

The new count-enhanced cyclic certificate is independently audited over full domains. Its published near-neighbor count inputs remain external mathematical dependencies, not newly replayed spherical SDP certificates.

## Guarded five-dimensional refinement

The current strict 5D upper is **0.510129239**, versus the same paper-reported 0.5122645 and D5 lower pi^2/(15*sqrt(2)). It lowers the pinned upper by **0.416828%** and closes **4.542443%** of the additive gap. Earlier 5D proofs and their original outward bounds remain preserved.

Only the guarded exact witness is certified. Its unguarded numerical checkpoint fails a far-domain condition and is not a valid bound. The exact positive-semidefinite kernel correction is explained and checked in the package; this is not a Gaussian sign repair. Published near-neighbor count inputs remain external dependencies.

## Five-dimensional lens extension

The current strict upper is **0.509996147**, compared with the same paper-reported 0.5122645 and D5 lower pi^2/(15*sqrt(2)). This reduces the pinned upper by **0.442809%**, and the additive gap by **4.825576%**. Prior 5D witnesses remain available and replayed.

This witness genuinely violates the former pointwise far-sign condition. A new long-edge lens budget is therefore essential: elementary centroid caps bound common neighbors, a verified physical-tail allowance pays the positive far contribution, and exact pair/triangle incidence accounting preserves the PSD density inequality. The independent audit treats step endpoints and periodic translation orbits explicitly. No new external spherical-code bound is used for the elementary lens caps; the existing near-neighbor count theorems remain external dependencies and their original SDP certificates are not reconstructed. No exhaustive literature-priority claim is made.

## Stronger five-dimensional spherical-lens follow-on

The current strict 5D upper is **0.509953498**. Relative to the unchanged paper-reported 0.5122645 and D5 lower, it decreases the upper by **0.451135%** and closes **4.916305%** of the additive gap. The earlier centroid-only lens bound 0.509996147 remains preserved.

This version strengthens the internal lens caps through conservative transverse projection and exact spherical LP certificates, combined with Rankin and centroid arguments. Those new cap certificates are included and replayed; they are not additional unverified external spherical-code inputs. The previously cited near-neighbor count theorems remain external dependencies. Both full primary and separate independent proofs check the changed physical-tail budgets, all boundary conventions and the relaxed far domain. The comparison remains pinned rather than an exhaustive current-worldwide-priority claim.

## New 43D cyclic bound

Strict upper 0.000002629078004 is compared with pinned 0.0000027169074727 (Cohn indexed table, citing 2006.02560v2). The newer lower from 2607.20359v4 is pi^21*137^21*sqrt137*22!/(2^67*sqrt2*3^21*44!), approximately 4.6423267938841084e-8. Upper reduction 3.232700%; additive-gap reduction 3.288897%. The lower is the new six-translate antipode construction, not the older table lower. No exhaustive current-worldwide priority claim is made.

## New 42D cyclic bound

Strict upper 0.000003677844262 is compared with pinned 0.0000038780970907. The newer lower from arXiv:2607.20359v4 is pi^21*3^17*sqrt3/(2^21*21!), approximately 5.7516899055409684e-8. Upper reduction 5.163688%; additive-gap reduction 5.241425%. Both comparisons use full decimals, not displayed approximations.

## New 41D cyclic bound

Strict upper 0.000005243057366 is compared with pinned 0.0000055306464395. The newer lower from arXiv:2607.20359v4 is 2^21*pi^20*3^17*21!/42!, approximately 8.6368303579508264e-8. Upper reduction 5.199918%; additive-gap reduction 5.282410%. Both comparisons use full decimals, not displayed approximations.

## Certified midrange continuation

The new25,28,32,34,36,38,39,40D entries use the pinned Cohn indexed upper table citing arXiv:2006.02560v2. Lower constructions are cited individually below; these are published inputs, not independently reconstructed packing constructions. Gap percentages use the full formulas and outward upper thresholds.

- 25D: L = 2^26*pi^12*13!/(26!*sqrt(2)); source https://cohn.mit.edu/sphere-packing/. New strict upper 0.001375214329.
- 28D: L = pi^14/14!; source https://cohn.mit.edu/sphere-packing/. New strict upper 0.0004827696028.
- 32D: L = pi^16*3^16/(16!*2^24); source https://cohn.mit.edu/sphere-packing/. New strict upper 0.0001249158969.
- 34D: L = pi^17*3^16*sqrt(3)/(17!*2^25); source https://cohn.mit.edu/sphere-packing/. New strict upper 0.00006394047941.
- 36D: L = 2^18*pi^18/(3^10*18!); source https://doi.org/10.1109/18.119683. New strict upper 0.00003121945691.
- 38D: L = pi^19*3^15*sqrt(3)/(19!*2^19*5); source https://arxiv.org/html/2607.20359v4. New strict upper 0.00001511771069.
- 39D: L = 2^20*pi^19*3^16*20!/(sqrt(5)*40!); source https://arxiv.org/html/2607.20359v4. New strict upper 0.00001064953048.
- 40D: L = pi^20*3^16/(2^20*20!); source https://arxiv.org/html/2607.20359v4. New strict upper 0.000007516663425.

## Gaussian cyclic continuation in 22, 23, 26 and 27 dimensions

These four entries improve the pinned Cohn upper table citing arXiv:2006.02560v2. The known lower constructions below are cited from https://cohn.mit.edu/sphere-packing/; the indexed snapshot was four months old and direct retrieval returned 502. They are comparison inputs, not independently reconstructed packing constructions or exhaustive current-record claims.

- 22D: known lower pi^(22/2)/Gamma(22/2+1) * (2^(-23) * 3^(-21/2) * 11^11); strict new upper 0.003653383161.
- 23D: known lower pi^(23/2)/Gamma(23/2+1) * (2^(-34) * 3^(-12) * 23^(23/2)); strict new upper 0.002642540089.
- 26D: known lower pi^13/(13!*sqrt(3)); strict new upper 0.0009746292767.
- 27D: known lower 2^28*pi^13*14!/(28!*sqrt(2)); strict new upper 0.0006864321082.

## Remaining 25–40D Gaussian cyclic certificates

The 29,30,31,33,35,37D upper comparisons are pinned Cohn-table values associated with arXiv:2006.02560v2. The following lower formulas are published construction inputs from https://cohn.mit.edu/sphere-packing/, not independently rebuilt packings. Individual independent audits record source review and scope.

- 29D: lower 2^30*pi^14*15!/(30!*sqrt(2)); strict new upper 0.0003420173821.
- 30D: lower pi^15/15!; strict new upper 0.0002454113168.
- 31D: lower 2^8*sqrt(2)*pi^15*3^15*16!/32!; strict new upper 0.0001763894399.
- 33D: lower 2^9*pi^16*3^16*sqrt(3)*17!/34!; strict new upper 0.00009009682341.
- 35D: lower 2^37*sqrt(2)*pi^17*18!/36!; strict new upper 0.00004463661144.
- 37D: lower 2^40*sqrt(2)*pi^18*19!/38!; strict new upper 0.00002171160843.

## 21D Gaussian refinement

Strict new upper0.004985977000 improves the prior repository upper0.005149549312 by3.176439% and closes7.369051% of its gap to the same known lower. Against the pinned published0.0051596603948176, the reductions are3.366179% and7.789082%. The existing lower formula/source is unchanged. See the [independent comparison](verifications/dim21/refinement_gaussian/independent/comparison.json); source-freshness caveats remain.

## Rebuilt angular kernels in 46D and47D

The radial data are unchanged; new target-dimensional angular kernels reduce the origin cost. Known lower formulas and pinned published upper references remain unchanged. Exact new upper bounds and comparison arithmetic are recorded in [release comparisons](research/angular_border_release/COMPARISONS.json). These are pinned-reference improvements, not exhaustive current-priority claims.

## Joint angular/radial refits in 46D and47D

The new strict bounds are0.0000008868492132 and0.0000006010083731. They improve the pinned published uppers by4.590129% and7.405509%, closing4.735895% and7.673514% of the corresponding additive gaps. Lower formulas and citations are unchanged. [Release comparisons](research/angular_radial_refit_release/COMPARISONS.json) also distinguish gains over the immediately preceding angular-only certificates.

## 15D curved-profile comparison

The pinned newer upper .03372585 is the paper-reported15D entry of [2512.04829v2, Table 1](https://arxiv.org/html/2512.04829v2#S3.T1), checked 2026-11-11; its external certificate is not replayed here. The known lower uses center density2^(-9/2), from the [Cohn table](https://cohn.mit.edu/sphere-packing/) indexed four months earlier, giving ordinary density L=2^12*pi^7*8!/(sqrt(2)*16!)=.0168575706567626... . Multiplication by unit-ball volume is recomputed independently at 70 and100 decimal digits. This is a pinned comparison, not exhaustive worldwide priority verification.
