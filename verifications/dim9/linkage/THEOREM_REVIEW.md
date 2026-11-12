# Independent analytic linkage review

Reviewed the frozen live2 construction; this review does not replace independent exhaustive geometry or radial sign checks.

The signed hinge condition is sufficient. On each interval g(v)=v phi'(v)=Av+B; A>=0, nonnegative jumps and g(1)>=0 imply g>=0 and nondecreasing. Thus phi is continuous, nonnegative and locally Lipschitz. For z>=y>=1, phi(xz)-phi(xy) increases with x, because its derivative is (g(xz)-g(xy))/x>=0.

For a sorted triangle with squared side ratios1,y,z, the three angular arguments are ordered t<=c1<=c2. Q(t)<0, and Q has exactly one negative-to-positive transition on(1/2,1). After dividing the cyclic sum by h(x)h(xy), its other multipliers are r1=h(xz)/h(xy) and r2=h(xz)/h(x). If both Q(c1),Q(c2) are nonnegative they decrease. If both are negative the conclusion is immediate. In the mixed case factor r1[Q(c1)+(h(xy)/h(x))Q(c2)]. The bracket decreases; while positive its product with decreasing positive r1 decreases, and once nonpositive it cannot defeat Q(t)<0. This justifies reduction to the compact/tail test at smallest scale.

For u>=4, y>=16. For t>=-1/4 the middle angle has cosine<=1/2 and is negative, leaving at most B exp(-phi(16)). For t<-1/4, z-y=1-2ut>=3. After all knots (<=2), phi'(v)=A+B0/v is >=gamma=A+min(B0,0)/16 for v>=16. The other positive contributions are bounded by B[exp(-3gamma)+exp(-phi(16))]. The exact structure check must bound this by epsilon.

The rooted kernel degree-zero block is C times(1,-1/B)(1,-1/B)^T, and all other Gegenbauer blocks are PSD. Pair contribution is -C/B(2h-h²). Cyclic contributions are nonpositive. Adding a pair slack lambda times the half-open near-neighbor indicator gives rho<=P(0)+C+373lambda, using average ORDERED neighbors per center, not a pointwise count, and without integer flooring.

I directly checked arXiv:1206.2608 Table2 contains dimension9 (epsilon,U)=(0.003013,373). The chosen threshold1.003012 is strictly inside it. The original published SDP is not replayed and remains an explicit external theorem dependency.

The Gaussian polynomial is Schwartz. The profile has positive final linear slope and logarithmic correction, so rooted sums converge absolutely. Standard periodic approximation and Poisson summation are analytical inputs rather than code-formalized results.

Independent normalization check uses min squared separation9/8, sphere radius squared9/32, unit9-ball volume32*pi^4/945. It verifies the coefficient linkage and objective enclosure below0.190926892, not the geometry or Fourier inequalities themselves. See verify_linkage.py and LINKAGE.json.
