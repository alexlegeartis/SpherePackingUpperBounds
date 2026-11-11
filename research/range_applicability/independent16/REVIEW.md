# Independent diagnostic16D mixed-Gaussian audit

This is a valid-but-noncompetitive applicability test, not a new record. The target saturated upper0.025150156 exceeds the newer published comparison0.02492121 reported in arXiv2512.04829v2.

The frozen source has exact exclusion square2 and truncation square14/5, seven radial hats, angular degree16 and Gaussian polynomial degree31. The source reciprocal-guard schema is checked to have zeta=0 and no Fourier repair intervals. It therefore reduces exactly to the independently reviewed slow-Gaussian lemma, with beta1/4, c157/25 and eta53/10^16.

In variables X=pi|x|²,Y=pi|xi|² the added physical term eta(c−X)exp(−X/4) has exact Fourier transform eta 2^16(c−32+16Y)exp(−4Y). Its negative Fourier part occurs only below Y=643/400. Exact Sturm tests prove H−eta 2^16(32−c)>0 there, while H is positive everywhere and the correction is nonnegative afterward. The corrected Fourier origin is exactly499999995533197963139/500000000000000000000. The physical origin gains eta c. The physical correction is negative at every allowed separation because2*pi>c.

The independent implementation checks exact rational Bernstein upper bounds on both physical repair intervals [59247/1000,59299/1000] and[14849/100,149581/1000], with exact Sturm negativity on every complementary tail interval. The first repair has a relatively small but strictly positive margin: polynomial upper about5,549,953.396 is below allowance5,575,259.960. These bounds use exact rational coefficients and Arb enclosures, not numerical sampling.

The original kernel pair inequality passes26 independently bounded leaves, maximum depth6. Exact rational elimination checks120 positive pivots. The direct monomial Fourier identity uses parameter n/2−1=7. All56 compact cyclic and21 far domains are checked through explicit solid harmonics and rational Bernstein subdivision. Infeasible far boxes may be removed only through the exact triangle-inequality test w>u+v and(w−u−v)^2>4uv at conservative rational box endpoints. The standard periodic cyclic proof is unchanged.

The unsaturated density enclosure is centered at0.0251948365559348121937635977395565217523462210715159078310748273632324436927925139859343975779801788301021301512678. Independent exact origin-block Schur factorization and rational Machin rounding prove saturated density below0.025150156, with enclosure centered at0.02515015598684426001957270545311974925532451731023909968650371208416686196761174108620928677618171167954130162714993.

Full verification passed all77 domains, using5 exact rational infeasible-triangle exclusions. See full.log and domains_audit.json. This test shows the same mechanism applies below17D; it does not show improvement in16D or universal optimizer success.
