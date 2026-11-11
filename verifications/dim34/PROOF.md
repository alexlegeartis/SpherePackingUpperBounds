# Proof linkage for the 34D certificate

## Claim and normalization

Let s=17/4, kappa=3, alpha=kappa*pi*s=51*pi/4>40. The two rational degree-111 polynomials P,H satisfy exactly

 f(x)=P(pi*|x|^2) exp(-pi*|x|^2),
 fhat(y)=H(pi*|y|^2) exp(-pi*|y|^2)

in physical dimension 34. The replay reconstructs this Fourier identity, proves H>0 on the complete nonnegative ray and P<0 on [7*pi,infinity), and checks H(0)>0. Root isolation covers all111 roots with multiplicities.

Let Q be the rational dimension 34 angular polynomial. Its normalized34D Gegenbauer coefficients are all nonnegative, its constant coefficient is1, and B=Q(1). Put h(r)=exp[-kappa*pi*(r^2-s)] and use the positive rational M in certificate.json:

 K(x,y)=(M/B) h(|x|)h(|y|)Q(angle(x,y)),
 K(x,0)=-M h(|x|), K(0,0)=MB.

The degree-zero origin Gram is rankone since (MB)(M/B)=M^2. All higher harmonics are PSD by the dimension 34 addition theorem. The pair condition f(r)+2K(r,0)+K(r,r)<=0 is rigorously checked on [sqrt(s),sqrt(7)]; beyond it follows from P<0 and0<h<=1. A smooth cutoff below sqrt(s) defines rapidly decreasing kernel features globally without changing packing-domain values.

The density conclusion is

 Delta_34 <= pi^(34/2)*(s/4)^(34/2)/Gamma(34/2+1) * (P(0)+MB)/H(0)
              < 0.00006394047941.

A rational Machin enclosure for pi proves the strict displayed decimal without relying on printed floating-point digits. The pinned comparator is6.45221967438e-05 from Henry Cohn's sphere-packing table, associated with arXiv:2006.02560v2. This package does not establish worldwide priority or optimality.

## Complete triangle domain

The auxiliary angular verifier proves Q<-1e-5 on[-1,1/2], Q<0 through5001/10000, and exactly one subsequent root before1. Write a triangle with shortest side1, second length u>=1, and included cosine t<=1/2, with2tu<=1. The third side is v=sqrt(1+u^2-2tu). Its cyclic sum after removing a positive Gaussian factor is

 Q(t)+exp[-alpha(1-2tu)]Q((1-ut)/v)
       +exp[-alpha(u^2-2tu)]Q((u-t)/v).

The complete compact region0<=t<=1/2,1<=u<=5 is covered by rational dyadic rectangles and Arb enclosures. Only rectangles satisfying2tu>1 throughout are excluded. The single-sign-change structure extends the proved alpha40 inequality to every alpha>=40 and every larger shortest side. The regions t<=0 and u>=5 are covered by the exact exponential tail inequality40^30/30!>2B*100000. Degenerate triangle boundaries are included.

This is the cyclic finite-sum refinement of the Cohn–de Laat–Salmon framework (arXiv:2206.15373), not the individual pointwise-sign premise of Theorem1.1. Absolute convergence permits grouping into cyclic triangle sums for periodic packings; periodic approximation gives unrestricted packing density. No lattice restriction or unproved neighbor-count assumption is used.

The exact angular and complete-domain replay scripts are frozen in continuum_three_point/. Counts and the final rigorous interval are reported by verify.py.
