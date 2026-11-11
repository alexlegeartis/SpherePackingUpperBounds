# Independent review: mixed-Gaussian wide-annulus 17D certificate

Frozen source: certificate17_repaired.json from research/density_solid_cyclic_20261111. SHA256 1bbc96f096ad55b58f2b63276618aa1350798589d09e0e6cf26687cef1900405. This review retains the exact exclusion square 17/8 and truncation square 29749999999999996/10^16, not the rounded value 119/40. There are seven radial hats and angular degrees 0–16.

## Independent mixed Fourier identity

Write X=pi|x|²,Y=pi|xi|². The physical correction is eta(c−X)exp(−beta X), with beta=1/4,c=667/100,eta=13/10^15. Differentiating the Gaussian Fourier transform with respect to beta gives its Fourier transform

eta beta^(−17/2)[c−17/(2 beta)+Y/beta²]exp(−Y/beta)
= eta 2^17(c−34+16Y)exp(−4Y).

Thus the corrected Fourier origin is H(0)−eta 2^17(34−c), checked exactly as 499999976715716832939/500000000000000000000. The physical origin gains eta c. Both normalizations must be included in the objective. The result is still Schwartz, so periodic Poisson summation remains available.

Below Y=2733/1600 the Fourier correction divided by exp(−Y) is at least minus cost=eta 2^17(34−c). Exact Sturm sequences verify H−cost>0 on this interval. Beyond that cutoff the Fourier correction is nonnegative, and the original H is globally positive. This proves the full mixed Fourier transform is positive.

The correction is strictly negative for every allowed physical separation because pi*(17/8)>667/100. Consequently the original verified annular pair inequality remains sufficient. It does not need to be reoptimized or checked by sampling.

## Full physical tail repair

The independent verifier proves P<0 between and outside the two rational repair intervals [2113/50,423/10] and [6109/100,1223/20] using exact Sturm sequences and negative endpoint values. It starts at the conservative rational cutoff R²*pi_lower. Within each repair interval it computes exact rational Bernstein coefficients for P. Their maximum is compared, using Arb, with the lower bound eta(lo−c)exp((1−beta)lo). This proves P(X)+eta(c−X)exp((1−beta)X)<0 over both intervals. These checks independently reconstruct the bounds rather than using the producer's Taylor enclosure routine.

## Cyclic argument, domains and PSD

The independent algebra verifier checks 120 positive rational pivots, the direct polynomial-Gaussian Fourier identity, and the original coupled annular pair inequality with 29 leaves of depth at most 6. Independent explicit rising-factorial solid harmonics and fully rational Bernstein subdivision cover 56 compact cyclic and 21 far cells.

For a far box, geometry exclusion is allowed only when exact rational endpoints satisfy w_min>u_max+v_max and (w_min−u_max−v_max)^2>4u_max*v_max. This proves w_min>(sqrt(u_max)+sqrt(v_max))² and excludes only triangles impossible by the triangle inequality. Boundary equality is never excluded. No floating-point value certifies a sign or geometry exclusion; floats only choose which coordinate to subdivide.

Use the same explicit periodic cyclic proof as the narrower 17D certificate: PSD sums around each center expand into the origin term, coupled pair terms, and twice unordered cyclic triangle sums. One long edge leaves one far kernel term, while two long edges leave none. The compact and far checks therefore make the total triple contribution nonpositive. Poisson summation yields number density at most the corrected objective; multiply by ball volume at radius sqrt(17/8)/2. Periodic approximation extends to arbitrary packings. This is not the literal pointwise-sign hypothesis of Theorem 1.1 and imposes no new packing assumption.

## Exact Schur saturation

For the origin block A=[[a,b^T],[b,B]], solve Bx=b exactly and set alpha=a−b^Tx. The verifier checks alpha>0 and reconstructs A−alpha e0e0^T=W B W^T with W=[x^T;I]. Hence saturation preserves PSD, all non-origin kernel values, the pair inequalities and triangle constraints. Only K(0,0) changes, so the corrected objective decreases by alpha/corrected_normalization.

Exact rational Machin bounds prove the unsaturated bound below 0.018250915 and saturated bound below 0.018204096. The saturated Arb enclosure is centered at 0.01820409584345106969660176457955375800073008402490479718223880953900950715779919142174730669776257801285864963631007. The full independent domain run passed all 56 compact cyclic and 21 far cells. Exactly three exclusions were used, each proved impossible by the rational triangle test. See full_exact_geometry.log and domains_audit.json.

This is a reproducible computer-assisted review, not a proof-assistant formalization or external peer review. The repair and Schur arguments change the certificate, not the underlying cyclic density theorem. No worldwide-priority claim is made.
