# Independent46D radial perturbation audit

The source is the frozen radial_perturb46_degree15_certificate.json copied to certificate.json. It contains the final exact degree111 physical and Fourier monomial polynomials, positive M and angular constant b0. Optimizer output, the degree15 discovery representation and rejected candidates are not proof inputs.

##46D radial pair

Unlike the47D restriction construction, this changes the radial function. The independent verifier computes all112 Fourier coefficients directly in dimension46 using the normalized Gaussian monomial/Laguerre formula with parameter22 and checks equality to the frozen H. The positive Fourier origin is H(0), not a48D normalization or a numerical optimizer estimate.

Exact rational Sturm sequences separately establish H>0 on[0,infinity) and P<0 on[7*pi_lower,infinity), with negative/positive endpoint signs as applicable. This is a different global-sign implementation from the producer's certified complex-root isolation. No small sign violation is tolerated, and no Gaussian correction is assumed.

Let s=pi|x|², k=17/8 and C=51*pi/4. The required raw coupling is P(s)exp(−s)<=M(2h−h²), h=exp(−ks+C). Dividing by positive h gives
F(s)=P(s)exp((9/8)s−C)−2M+M exp(−ks+C)<=0.
The independent verifier encloses this complete function on the slightly enlarged rational interval[6*pi_lower,7*pi_upper] using F(mid)+|F'(mid)|radius+sup|F''|radius²/2. The polynomial part of the second derivative is independently reconstructed as P''+(9/4)P'+(81/64)P. Centered Arb polynomial enclosures and monotonic exponential bounds prove each subdivision. Beyond physical radius squared7, P<0 and0<h<=1, so coupling follows immediately.

## Kernel and cyclic dependencies

The exact degree19 angular source Q is unchanged, with B=Q(1). Explicit rising-factorial polynomials reconstruct its46D normalized Gegenbauer expansion, check every coefficient nonnegative, and separately verify b0 by spherical moments. At allowed nonzero arguments the kernel remains(M/B)h(x)h(y)Q(angle). The origin/degree-zero Gram matrix has entries MB/b0,−M,M*b0/B and determinant0. A smooth radial cutoff gives the same smooth PSD extension as in the audited47D construction. Thus K(x,0)=−Mh, K(x,x)=Mh², and K(0,0)=MB/b0.

All non-origin cyclic values depend only on Q, k and the three side lengths. They inherit the audited48D global cyclic sign, including compact rectangles, analytic tails and the scale lemma. Fresh source angular and domain scripts are replayed in this directory; their logs are required proof dependencies. Minimum separation squared remains6, and(51/4)*pi>40 is independently checked. No lattice-only or additional packing assumption is used.

The periodic PSD/cyclic argument and Poisson summation yield density at most the46D ball volume at radius sqrt6/2 times(P(0)+MB/b0)/H(0). The volume is(3*pi/2)^23/23!. A rational Machin upper bound certifies the final outward decimal0.000000914882847, separately from Arb display arithmetic.

## Status and comparison

All new radial gates and both inherited source gates passed. Exact Sturm variations are55 to55 for both Fourier and physical tails. The independent coupling bound passes10 leaves with maximum depth8. The final density enclosure is centered at9.1488284600262300906264110774671012774815093433564050608912408891658744426658146413084020013286703491110482354621588645003319793508420592842195065235930e-7. See audit.log/audit.json and both source replay logs.

Cohn's indexed46D table row gives pinned upper0.0000009295151556. The available index was four months old and the live site returned502 during this audit. The comparison is not a guarantee of exhaustive current worldwide priority. This is a reproducible computer-assisted audit, not external peer review or proof-assistant formalization.

Frozen certificate SHA256: d26fb1e751f1a475d1955531ef04b8e120a1827299707bb37fd1f4ae1aeb1aaa.
