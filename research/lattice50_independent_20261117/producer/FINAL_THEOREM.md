# 50D lattice-only upper certificate

Producer certificate passed; independent complete audit is pending. This is not a comparison against a verified literature record.

Use the exact rational Fourier-positive polynomial Q_F and spatial polynomial P in repaired.json, with f(x)=exp(-x)P(x), Fourier convention x=pi|v|² and minimum x22. Q_F(0)=1. Positive rational linear/quadratic factors establish Q_F>=0 on[0,infinity), and exact Laguerre inversion proves the Fourier identity.

Use the earlier independently audited scalar rank-two angular kernel, positively rescaled to constant Gegenbauer coefficient1 in dimension50. All exact50 Gegenbauer coefficients are nonnegative; its contact q=1150750294032540977074575/1437673731494371. Its h(t)=t^96 exp[-120(t-1)] rank-two sign theorem is dimension-independent and unchanged by this positive scalar rescaling. The old48 radial witness is not used.

Let a and z be the exact rationals in annular_final_certificate.json and set s=z²/a+10^-12. Thus [[s,z],[z,a]] is positive definite. Let rho=10^-70 and L=(1+rho)/(1-rho). The primitive-line inequality (sum h_k)^2 <= L sum h_k² holds from the audited geometric ratio bound. The pointwise corrected function is
 R(x)=exp(-x)P(x)+(2/q)[a L h(x/22)^2+z h(x/22)].

The final verifier proves R(x)<=g_i on every original annulus[x_i,x_i+1/4], from22 to180, using a combined Taylor polynomial and rigorous exponential remainder at768bits with adaptive subdivision. The recorded g_i are exact nonnegative rationals. The annular coset lemma in ANNULAR_CAP_LEMMA.md gives at most C_i unoriented pairs per mod-two coset in each annulus. All2^50 cosets, including zero, are included. Therefore sum over nonzero lattice vectors of g is at most B=2^51 sum C_i g_i.

Beyondx180 the exact physical polynomial is negative by the separately verified complete root isolation and sign test. Also h<=1 and aL+z<0, so the kernel correction is nonpositive there. Thus g=0 beyond180.

Construct K(0,0)=s, K(0,v)=z h(v²)/q, K(v,w)=a Q_A(cos(v,w))h(v²)h(w²)/q² for nonzero vectors. Its total lattice sum is nonnegative by PSD. The independent-pair sum is nonpositive by the previously audited rank-two theorem. Combining the remaining primitive-line terms with R<=g yields sum_{v!=0}f(v)<=s+B. Poisson summation gives the packing density bound
 (22/4)^25 /25! *[P(0)+s+B].

All functions are Gaussian times polynomials, so absolute convergence and Poisson summation hold. Exact rational normalization in density_result.json gives strict safe upper1.75102693e-7. This number includes the10^-12 Schur guard and the actual exact annular heights; the earlier normalization_check.py's smaller10^-20 guard is not used for final density.
