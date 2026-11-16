# Cyclic PSD kernel plus lattice profile budgets: independent bridge audit

Normalize lattice minimum to at least 1. Let K(x,y) be a real symmetric Schwartz positive-semidefinite kernel on R^n x R^n. Assume, for all triangles with x,y,x-y nonzero and lengths at least 1,

 K(x,y)+K(-x,y-x)+K(-y,x-y)<=0.

Let C=K(0,0), D(x)=K(x,0)+K(0,x)+K(x,x). Let f be real Schwartz with nonnegative Fourier transform and positive integral fhat(0). Suppose for every |x|>=1,

 f(x)+D(x)<=G(x),

where G is radial, nonnegative, compactly supported below squared radius 4. Suppose the independently certified lattice class-budget theorem gives sum_{x in Lambda, x!=0} G(x)<=2(2^n-1)kappa, using the chosen mod-2 short/long profile constraints.

Then

 rho(Lambda) <= [f(0)+C+2(2^n-1)kappa]/fhat(0).

## Proof and multiplicities

The full kernel sum S=sum_{x,y in Lambda}K(x,y) is nonnegative: each finite square sum is PSD, and Schwartz decay gives absolute convergence to S. Split

 S=C+sum_{x!=0}D(x)+T,
 T=sum_{x,y!=0,x!=y}K(x,y).

Each of the three cyclic transformations is a bijection of the off-diagonal ordered lattice-pair domain. Thus 3T is the sum of the cyclic expressions above and T<=0. Consequently C+sum D=S-T>=0. This is an exact infinite lattice identity; no periodic averaging or boundary error is required.

Poisson summation yields rho*fhat(0)<=sum_{x in Lambda}f(x). The pointwise inequality and sum D>=-C yield

 sum f<=f(0)+sum_{x!=0}G(x)-sum_{x!=0}D(x)<=f(0)+C+sum G.

Combine the inequalities. The origin contributes C exactly once, not twice. The factor 2 in the profile term comes solely from the two signs in each antipodal pair. Cyclic triples are already accounted for by 3T and introduce no extra factor.

## Published radial-harmonic kernel specialization

With normalized Gegenbauer expansion Q=sum c_l P_l, c_0=1, c_l>=0 and B=Q(1), the allowed-domain kernel is

 K(x,y)=C/B^2 h(|x|^2)h(|y|^2)Q(cos angle),
 K(0,x)=-C/B h(|x|^2), K(0,0)=C.

The origin plus constant-harmonic coefficient matrix [[C,-C/B],[-C/B,C/B^2]] is PSD for C>=0. Other harmonic coefficients are nonnegative. A smooth radial cutoff below radius 1 gives a global Schwartz feature realization without changing any allowed-domain value. Hence

 D(x)=C/B [h(|x|^2)^2-2h(|x|^2)].

For h(v)=v^96 exp[-120(v-1)], 0<h(v)<=1 on v>=1. Therefore D<=0, and f<=0 beyond a certified cutoff suffices for the combined pair tail condition with G=0. The cyclic inequality remains a separate prerequisite supplied by the published angular/profile certificate; the displayed diagonal formula alone does not prove it.

## Optimization and normalization precautions

C must be constrained nonnegative. The origin value, Fourier convention, and any rescaling must remain consistent. If f_pub is defined on a packing of minimum sqrt(s), the min-one function is f_tilde(x)=f_pub(sqrt(s)*x), with fhat_tilde(0)=s^(-n/2)fhat_pub(0). The density is rho times volume of the radius-1/2 ball.

For an actual baseline guarantee, include the published f_pub as an exactly representable column or verify its polynomial lies in the chosen radial basis; merely importing its kernel into a different finite radial basis does not guarantee reproducing the published bound. Taking G=0 then recovers that bound.

Status: mathematical bridge proved under the explicitly stated hypotheses. No new numerical density improvement asserted here. Full exact cyclic/angular source certificates, combined pair-condition enclosures, Fourier positivity and lattice profile budgets must all be replayed for a numerical claim.
