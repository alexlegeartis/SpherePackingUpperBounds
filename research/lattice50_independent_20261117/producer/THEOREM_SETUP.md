# Direct 50-dimensional scalar-kernel screen: conditional proof setup

This is not a certified bound. Global Fourier positivity and the radial inequality remain to be certified.

Let Q be the earlier certified scalar 48-dimensional rank-two angular polynomial, renormalized so its constant Gegenbauer coefficient in dimension50 is1. Its exact nonnegative dimension50 expansion and q=Q(1) are in ../scalar48_kernel_upward_extensions.json. Multiplication by a positive constant does not change its previously certified rank-two sign theorem with h(t)=t^96 exp(-120(t-1)). No new geometry gate is needed.

Take a>=0 and [[s,z],[z,a]] positive semidefinite. For a lattice with minimum norm at least1, use nonzero kernel K(x,y)=a Q(cos(x,y))h(|x|²)h(|y|²)/q², K(0,x)=z h(|x|²)/q, K(0,0)=s. Gegenbauer expansion and the origin block prove positive semidefiniteness. The sum over independent pairs is nonpositive by the earlier rank-two theorem.

Put rho=10^-70 and L=(1+rho)/(1-rho). The profile ratios along successive positive multiples are at most rho. Thus H² <= L sum h_k² for H=sum h_k. A sufficient pointwise radial condition is

f(r)+2a L h(r²)²/q+2z h(r²)/q <=0 (r>=1).

Summing twice over each primitive unoriented line yields 2sum f +4a H²/q+4zH/q<=0. Together with PSD and the rank-two sign, this gives sum_{x!=0}f(x)<=s. Poisson summation closes the density inequality.

The numerical solver uses x=pi R² with minimum x=d=22, and writes the Fourier-normalized scalar function as e^-x sum c_k L_k^24(2x)/binom(k+24,k). Its Fourier transform is e^-x sum (-1)^k c_k L_k^24(2x)/binom(k+24,k), with value1 at0. Rescaling minimum norm to1 multiplies the Fourier integral by (pi/d)^25. Hence the density objective is (d/4)^25/25! times (f(0)+s). The current solver omitted L-1 (~2e-70); exact repair MUST restore it.

Remaining requirements: exact positive Fourier polynomial, exact origin PSD, exact all-r>=1 radial sign with restored L, convergence, and independently checked normalization and prior geometry dependencies. A sampled solver objective is not a bound or a record.
