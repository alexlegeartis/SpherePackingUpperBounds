# Independently checked origin-ratio modification

Let Q have nonnegative normalized Gegenbauer coefficients and constant coefficient 1. Put B=Q(1)>0 and qminus=Q(-1). Let C>=0 and r>=0. For nonzero points use the unchanged kernel C h(x)h(y)Q(angle)/B^2. Set K(0,0)=r^2 C and K(0,x)=-r C h(x)/B. The constant-harmonic augmented block is

[[r^2 C, -r C/B],[-r C/B,C/B^2]],

which is positive semidefinite (an outer product). All remaining harmonic coefficients stay nonnegative. The cyclic inequalities involve nonzero triangle differences, so the already verified angular cyclic gate remains unchanged.

The lattice orbit sum over the three collinear ordered-pair types remains exactly the same. Consequently the per-vector correction is

D_eff(t)=C/B*((1+qminus/B)h(t)^2 - 2r h(t) + 2h(t)h(4t)).

Its lattice sum is at least -r^2 C, by the full PSD sum and cyclic grouping. Thus the density numerator is p(0)+r^2 C+2(2^n-1)kappa+extra. C itself is the nonorigin amplitude; replacing every C by r^2 C would be incorrect.

For the candidate r=7/10 is exact. The audit verifies origin_C=C*r^2, inserts r only into the negative cross term in the independent centered mixed-exponential enclosure, and uses origin_C in the objective. Since h decreases on t>=1, 0<h<=1, and a=1+qminus/B>0, the exact check a+2h(4)<2r proves D_eff<=0 throughout the allowed domain, including the polynomial tail. Gaussian convergence justifies all lattice sums.
