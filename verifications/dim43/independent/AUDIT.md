# Independent43D Gaussian cyclic certificate audit

PASS: unrestricted sphere-packing density <0.000002629078004.

Frozen certificate SHA256:1d8f2388a010d71ca71c605eba7fe1fc85acd8ee43603808202ce263ead3e821. The source-hash manifest binds both certificate.json and the new dimension43 angular polynomial. The standalone verify_all.py runs both independent proof implementations; python-flint is the only external verification dependency.

## Exact radial and origin checks

Minimum squared separation is45/8, physical-tail threshold7, kappa34/15, and kappa*pi*(45/8)=51pi/4>40. All112 monomial Fourier coefficients were independently reconstructed in dimension43 and match exactly. Exact rational Sturm sequences prove no roots of H on[0,infinity), with positive endpoint, and no roots of P on[7*pi_lower,infinity), with negative endpoint. Both variation counts are55→55. These checks use exact integer/rational polynomial arithmetic, independently of the producer's complex interval root isolation.

Explicit rising-factorial dimension43 Gegenbauer polynomials reconstruct Q. Every harmonic coefficient is nonnegative, their sum is B, and the constant coefficient b0=1 also agrees with an independent spherical-moment calculation. The nonorigin degree-zero coefficient M*b0/B and isolated-origin value M*B/b0 have product M², giving the exact rank-one Gram extension. Thus the origin cross term is-Mh, diagonal isMh², and the coupled pair condition is f<=M(2h-h²).

The independent coupled verifier uses second-derivative Taylor bounds on the complete divided expression. It proves all radii squared in[45/8,7] with17 closed leaves, maximum depth6. Beyond7, f<0 and0<h<=1 complete the ray. H(0)>0 and the exact objective (P(0)+MB)/H(0) match the source.

For dimension43 the radius sqrt(45/8)/2 ball volume is

2*pi^21*(45/8)^(43/2)*22!/44!.

Squaring it removes the square root:4*pi^42*(45/8)^43*(22!)²/(44!)². Exact rational upper Machin bounds prove objective² times this squared volume is strictly below(0.000002629078004)². The independent Arb density enclosure is centered at2.62907800321760827045287294452590036991782931052891848162494072593637004605022394261489200067569416295484986558796554278122650934991657842957942447920375e-6, radius below6.72e-159.

## Complete angular and cyclic checks

The independently implemented angular audit checks strict Q<-1e-5 on[-1,1/2], Q<0 on[-1,5001/10000], and exactly one sign-changing root thereafter before1. Its separate exact-Bernstein/Arb triangle implementation proves the complete sorted compact domain using5,554 rectangles, with421 exact sorting exclusions, maximum depth15. These exclusions have impossible sorting geometry, not omitted valid triangles.

Analytic negative-angle and large-side-ratio tails follow from the newly verified B and exact exponential-series bound. The single-sign-change scale lemma extends the alpha40 calculation to51pi/4 and all larger shortest edges. In the mixed-sign case the positive secondary term decays at least as fast as the negative one; both-negative and both-positive cases are immediate. The dimension-specific sign and coefficient hypotheses were checked on the new Q, not inherited from an older witness.

A smooth radial cutoff at the origin gives a globally smooth rapidly decreasing PSD kernel while preserving all required packing evaluations. Its cyclic triangle sum, rather than individual off-diagonal negativity, is the required hypothesis. Summed periodic PSD/Fourier bookkeeping and unrestricted periodic approximation are unchanged from the established cyclic framework. No near-neighbor count or new geometric external input is used.

## Comparison and defect corrected

The producer's initial proof note contained a wrong44D-like comparator .0000019017703143. It was reported and corrected to the pinned43D upper .0000027169074727 before integration. This was a documentation error; it did not enter the independent certificate checks.

The live primary paper https://arxiv.org/html/2607.20359v4, Table1 and Section4, supplies a newer43D lower construction than the older Cohn table index: center density137^(43/2)/(2^(223/2)*3^21). Therefore packing lower density is

L=pi^21*137^21*sqrt(137)*22!/(2^67*sqrt2*3^21*44!)
 =4.64232679388410840588823802085953247505…e-8.

Against pinned U=.0000027169074727 and outward N=.000002629078004, the upper decreases3.232700030550437%, and the additive gap to this newer lower closes3.288896767987259%. The lower construction is sourced from that paper, not independently reconstructed here. The comparison is not an exhaustive current-worldwide-priority claim.

No mathematical defect remains in the frozen certificate. This is independent computer-assisted checking, not proof-assistant formalization or external peer review.
