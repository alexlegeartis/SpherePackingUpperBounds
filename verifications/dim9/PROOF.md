# 9D signed-profile density bound: proof and ideas

The included exact primary checks and separate audits establish density <0.190941549, conditional on the explicitly cited published average-count theorem. This is a fresh reconstruction, not byte-identical recovery.
Use minimum squared separation s=9/8, v=r²/s and h(r)=exp(-phi(v)).
For the linear/log hinge profile, check exactly that g(v)=v phi'(v) is
nonnegative and nondecreasing: g(1)>=0, each cumulative linear slope>=0,
and every knot jump k*a_k+b_k>=0. Individual coefficients may be signed.
Then phi>=0 and h lies in(0,1]. For x>=1 and v>=1,
phi(xv)-phi(x) is nondecreasing in x because its derivative is
[g(xv)-g(x)]/x. Thus h(xv)/h(x) is nonincreasing.

Let Q have nonnegative normalized Gegenbauer coefficients with constant1,
B=Q(1). Check Q<=-epsilon on[-1,1/2] and exactly one transition from
negative to positive in(1/2,1). For sorted normalized triangle lengths
1<=u<=w, w²=1+u²-2ut, prove the signed expression
Q(t)+exp(phi(u²)-phi(w²))*Q((1-ut)/w)+exp(-phi(w²))*Q((u-t)/w)<=0.
Compact coverage is t∈[-1,1/2],u∈[1,4],2tu<=1. Full rectangle leaves
must partition exact area9/2, with excluded leaves wholly infeasible.
Signed interval mean-value evaluation remains valid across hinge knots
because the profile is locally Lipschitz and derivative intervals include
both one-sided values.

For u>=4 set final slope lower gamma=A+min(B0,0)/16, where final profile
is A*v+B0*log(v)+constant. Verify gamma>=7. Tail contributions are bounded
by B*exp(-phi(16)) if t>=-1/4, or
B*(exp(-3gamma)+exp(-phi(16))) otherwise; both must be<epsilon.
Scale reduction uses monotonic ratios. Both later angular terms negative
is automatic. If both are positive, both multipliers decrease. In the
mixed case their ordering forces middle negative and last positive;
factor the middle multiplier, and the positive last-to-middle ratio
also decreases. While the bracket is positive its product with the
decreasing positive middle multiplier decreases; once the bracket is
nonpositive, its contribution is nonpositive and Q(t)<0 finishes the case.
This covers all scales.

Rooted kernel K(0,0)=C, K(0,x)=-(C/B)h(|x|),
K(x,y)=(C/B²)h(|x|)h(|y|)Q(cos(x,y)) is PSD: its constant harmonic block
is C*(1,-1/B)^T*(1,-1/B), and all other blocks follow the addition theorem.
Cyclic triangle contributions are nonpositive by the previous check.

For f(x)=exp(-pi|x|²)P(pi|x|²), verify exact Fourier partner H using
normalized Laguerre polynomials L_k^(7/2)(2X)/L_k^(7/2)(0), signs(-1)^k,
H(0)=1, H>=0 everywhere and P<0 for X>=7. Verify every point of the pair
inequality f(x)<=(C/B)(2h-h²)+lambda*1[1<=v<(1.003012)²], lambda>=0.
Split at every hinge and count threshold, testing no-count side at threshold.
For v>=2, pi*s*v>7 and h∈(0,1] implies the pair automatically.

Published average ordered-neighbor count373 supplies at most373*lambda.
Poisson summation and summed PSD rooted forms imply center density
rho<=P(0)+C+373*lambda. Gaussian tails ensure convergent periodic sums;
standard periodic approximation extends to arbitrary packings. Multiply by
pi^(9/2)/Gamma(11/2)*(9/32)^(9/2), and certify outward rounding.
External analytic inputs: Cohn–de Laat–Salmon cyclic finite-sum framework,
https://arxiv.org/abs/2206.15373; addition theorem, Poisson summation,
periodic approximation, and the average-count theorem in EXTERNAL_COUNT.md.
The scripts are not a formalization of those analytic theorems.

## Certificate and audit results

The exact objective enclosure is 0.19094154845873147291964..., strictly
below0.190941549. Full compact primary geometry accepts20,362 boxes and
discards351, partitioning exact area9/2. Separate analytic-gradient replay
checks20,502 boxes after70 refinements of that partition. The degree63
radial proof uses51 intervals; a separate monomial Fourier/Sturm auditor
rechecks the full pair range using first-derivative bounds.

For signed profiles, the primary phi enclosure is obtained from endpoint
values, and phi' from [g(left−),g(right+)]/v, using exact monotonicity.
All normalized squared lengths are at least1 on the complete verification
rectangle, so clipping outward-rounded interval endpoints to1 is valid.
The separate geometry implementation uses natural signed expressions and
independently calculated analytic derivatives, with adaptive refinement.

Compared with the pinned paper-reported0.19099063 in arXiv:2512.04829v2
Table1, the bound is smaller. This does not independently certify that
paper's own rounding or establish worldwide priority or saturation.
