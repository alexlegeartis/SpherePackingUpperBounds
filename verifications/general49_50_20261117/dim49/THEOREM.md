# Unrestricted cyclic-kernel bound with annular error budget

This construction is for arbitrary sphere packings, NOT only lattices. It uses no parity classes, cosets, lattice pair-closure, or antipodal neighbor assumption.

Use f(v)=exp(-x)P(x), x=pi|v|², minimum x=22, and h(t)=t^192 exp[-216(t-1)], t=x/22. Exact Fourier polynomial Q_F is nonnegative on [0,infinity), Q_F(0)=1.

Let Q_A have nonnegative dimension-n normalized Gegenbauer coefficients, constant coefficient1, B=Q_A(1). Its complete cyclic triangle inequality is verified by the angular structure/infinite-tail and compact Taylor programs. For a>0 and s>=z²/a define
K(0,0)=s; K(0,v)=z h(t)/B;
K(v,w)=a Q_A(cos(v,w))h(t_v)h(t_w)/B².
The spherical addition formula and matrix [[s,z],[z,a]] give PSD. A smooth radial cutoff near the excluded origin realizes this as a global Schwartz PSD kernel, without changing allowed values. Its diagonal is a h(t)²/B, NOT twice that quantity.

Suppose R(x)=f(x)+(a h(x/22)²+2z h(x/22))/B <= g_i on annulus x in [l_i,r_i], with nonnegative exact g_i, and R<=0 beyond all bins. For any center in an arbitrary packing, directions of neighbors in this annulus form a spherical code with maximal inner product
u_i=1-22/(2r_i), provided (r_i-l_i)/22<=1.
Indeed their squared lengths divided by the minimum squared distance lie in [l_i/22,r_i/22]; pairwise distance>=1 gives cos<= (u+v-1)/(2sqrt(uv)), which increases in both variables when the bin width<=1.

For n>=2, any integer k>=0 with k/(k+n-1)>nu gives the explicit Delsarte polynomial
F(t)=(1+t)^k(t-nu), negative on [-1,nu]. Its constant Gegenbauer coefficient is
mu_k*(k/(k+n-1)-nu), mu_k=2^k*((n-1)/2)_k/(n-1)_k.
All other coefficients are nonnegative. To see this, the positive coefficients a_l(k) of (1+t)^k have ratio
 a_l(k+1)/a_l(k)=(k+1)(2k+n-1)/[(k+1-l)(k+n-1+l)],
which is at least the l=0 ratio 1+k/(k+n-1), since n>=2. The new highest-degree coefficient is positive. Therefore the number of neighbors is at most
C_i=floor[(1-nu)*(n-1)_k/(((n-1)/2)_k*(k/(k+n-1)-nu))].
There is NO factor2 and NO factor2^n in this unrestricted count.

Consequently every center has sum g<=E=sum C_i*g_i. Closed-bin endpoint overlap only overestimates this nonnegative sum. In the periodic PSD cyclic summation theorem, adding these per-center errors changes the bound from f(0)+K(0,0) to f(0)+K(0,0)+E. Explicitly the oriented distinct-triangle contributions are nonpositive by the cyclic condition, diagonal/origin terms contribute the displayed correction, and the PSD kernel sum is nonnegative. Poisson summation gives center intensity <=P(0)+s+E. Periodic approximation extends to arbitrary packing density.

Final density bound in dimension49 is V_49(sqrt(22/pi)/2) * [P(0)+s+E]. Its square is normalized exactly using a rational Machin bracket for pi in verify_density.py. This proof uses the cyclic version in Cohn-de Laat-Salmon Section2.2, equation(2.4), Theorem2.4; it does not substitute cyclic negativity into the stronger pointwise premise of Theorem1.1.

Certification status is recorded by separate outputs; this note alone is not a claim that every numerical check has completed.
