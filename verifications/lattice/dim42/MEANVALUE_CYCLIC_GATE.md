# Centered-gradient pruning for the complete cyclic domain

This adds a rigorous sufficient test to the existing outward-Arb rectangle verifier. It does not change the frozen angular polynomial or radial profile, the PSD check, exact negative cap, Sturm transition test, infinite tail, or all-dilation argument.

Write p=-a, b>0, z=1+u^2-2tu, w=sqrt(z), alpha=(1-ut)/w, beta=(u-t)/w, and

 R1=(z/u^2)^p exp[-b(z-u^2)],
 R2=z^p exp[-b(z-1)],
 F(t,u)=Q(t)+R1 Q(alpha)+R2 Q(beta).

The new test is used only if the entire rectangle lies in the sorting domain: max of its four endpoint products tu is <=1/2. The rectangle already has u>=1 and t<=1/2. Its center-to-point segments remain inside this rectangle and hence inside the sorting domain. Thus the existing monotone endpoint enclosures for R1 and R2, intersected with [0,1], are valid everywhere on every such segment.

The derivatives used are

 Lt = d_t log R1 = d_t log R2 = -2u(p/z-b),
 L1u = 2p[(u-t)/z-1/u]+2bt,
 L2u = 2(u-t)(p/z-b),
 alpha_t=-u/w+alpha*u/z,
 alpha_u=-t/w-alpha*(u-t)/z,
 beta_t=-1/w+beta*u/z,
 beta_u=1/w-beta*(u-t)/z.

Therefore

 F_t=Q'(t)+R1[Lt Q(alpha)+Q'(alpha)alpha_t]
                 +R2[Lt Q(beta)+Q'(beta)beta_t],
 F_u=R1[L1u Q(alpha)+Q'(alpha)alpha_u]
                 +R2[L2u Q(beta)+Q'(beta)beta_u].

Every expression is evaluated with outward Arb intervals over the full rectangle. Q' is reconstructed as the exact polynomial derivative; its interval enclosure uses the same centered Taylor composition as Q. The terms of each derivative are combined before taking an absolute-value upper bound.

For center (tc,uc), half-widths ht,hu, and interval derivative bounds Dt,Du, the multivariate mean-value theorem gives

 F(t,u) <= F(tc,uc)+ht sup|Dt|+hu sup|Du|.

Only when this outward upper bound is strictly negative is the entire rectangle accepted. Otherwise the previous verifier's split rule is followed. Partly invalid rectangles are never accepted by this new test; the old enclosure remains available on them. All accepted and discarded boxes still form the complete original subdivision tree.

All analytic derivatives and the whole-rectangle domain requirement were independently reviewed by the dual worker, with no missing chain or radial factors found. The first frozen41D eta=.15 angular kernel passed this augmented exact gate with23,498 accepted rectangles. This proves an angular/cyclic kernel only; a density claim still needs Fourier positivity, mixed radial envelopes, lattice profile and AA budgets, and normalization.
