# Independent45D new-angular and radial audit

This witness changes both the angular polynomial and radial Fourier pair. It is not the unchanged48D restriction construction. The frozen exact inputs are certificate.json (degree111 radial P,H,M,b0) and angular.json (new45D degree19 Q). Numerical optimizer output and rejected search candidates are not proof inputs.

## New angular proof

The independent algebra verifier reconstructs normalized45D Gegenbauer polynomials from explicit rising-factorial formulas with parameter43/2. Exact triangular expansion of Q has all coefficients nonnegative, constant coefficient b0=1, and sum B=Q(1). The spherical moment identity independently confirms b0. The addition theorem therefore proves PSD and |Q(t)|<=B on[-1,1].

verify_angular_cyclic.py independently proves Q<-10^-5 on[-1,1/2], Q<0 on[-1,5001/10000], and exactly one root between5001/10000 and1 using rational Sturm sequences. The positive endpoint Q(1)=B fixes the sign after that root.

For a triangle with sorted normalized sides1,u,v, put t=(1+u²−v²)/(2u), so t<=1/2,u>=1 and2ut<=1. The independent compact-domain enclosure covers t in[0,1/2],u in[1,5], discarding boxes only when2*t_min*u_min>1. It independently constructs exact Bernstein bounds for Q(t) and derivative-based Arb Taylor intervals for Q((1−ut)/v),Q((u−t)/v), v=sqrt(1+u²−2ut). It certifies Q(t)+Q((1−ut)/v)e^(−40(1−2ut))+Q((u−t)/v)e^(−40(u²−2ut))<0 on6,374 rectangles, with420 exact sorting exclusions and depth15.

The tails are mathematical, not sampled. For t<=0 both exponential gaps are at least1. For u>=5,t>=0, positivity of the first secondary Q forces its cosine>5001/10000, hence its gap>4.001; the second gap is at least24. The rational inequality40^30/30!>2B*10^5 bounds both positive tails by less than10^-5.

The scale lemma is unchanged but its premise is reverified for the newQ. The two secondary angles are ordered; a positive first term forces the second nonnegative. For exponents0<=p<=q, A*e^(-pT)+C*e^(-qT)<=max(A*e^(-p)+C*e^(-q),0) forT>=1 in the allowed sign patterns. Thus compact negativity extends to all shortest side lengths>=1 and all damping>=40. The actual damping51*pi/4>40 is independently checked.

## New radial proof and normalization

The direct45D monomial Fourier formula (Laguerre parameter43/2) checks all112 coefficients against H. Exact Sturm sequences, independent from producer complex-root isolation, prove H>0 on the nonnegative ray and P<0 beyond7*pi_lower.

The coupling function after division by positive h is F(s)=P(s)e^((9/8)s−51*pi/4)−2M+M e^(−17s/8+51*pi/4), s=pi|x|². An independent second-derivative Taylor enclosure proves F<=0 on the enlarged rational interval[6*pi_lower,7*pi_upper]. Beyond squared radius7, P<0 and0<h<=1 imply coupling. The exact positive normalization is H(0).

The kernel outside separation sqrt6 is(M/B)h(x)h(y)Q(angle). Its origin/cross degree-zero Gram entries are MB/b0,−M,M*b0/B, with positive diagonal entries and determinant0. A smooth radial cutoff supplies a smooth global PSD extension. Thus the required pair and cyclic inequalities yield the unrestricted periodic density bound, and periodic approximation covers arbitrary packings. No lattice-only or additional packing assumption is used.

The objective is(P(0)+MB/b0)/H(0). For radius sqrt6/2, the45D ball volume squared is4*pi^44*6^45*(23!)²/(46!)². Rational Machin pi bounds prove the outward target0.000001311423575 independently of Arb decimal display.

## Status and comparator

Both independent proof paths passed. Global Fourier/physical Sturm variations are55 to55 each. The new radial coupling passes12 leaves of maximum depth8. The density enclosure is centered at1.31142357458920838673180174230383870790282477601649408459057062169615871317630153434557563753711344822992790154575779308003074412056047300807789439711439e-6. See audit.json/audit.log and angular_cyclic_audit.json/angular_cyclic.log. The pinned Cohn45D comparator is0.0000013300905665; the indexed table was four months old. This is not an exhaustive worldwide-priority guarantee, external peer review, or proof-assistant formalization.

Fresh targeted literature search found no stronger finite45D general upper. arXiv2607.20359v4 improves the45D kissing number but explicitly ties the known density; it is not a better density comparator. The known center density remains23^(45/2)/2^(183/2), from Chen et al., arXiv2505.02394v1.

certificate.json SHA256: 5b62e11e7db1e2c5c94ce6f61ef9d0617fad6bbb62183534016845713661b757.

angular.json SHA256: 4ad86888f870b356df79112057fbffe136d6d04d9fb74af55beb9f476a79f750.
