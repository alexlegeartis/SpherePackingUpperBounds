# Independent analytic audit of live2 signed-profile construction

The signed profile is continuous and locally absolutely continuous. Its g(v)=v phi'(v) has affine slope A_i between knots, and jump k*b_k+d_k. Nonnegative interval slopes and jumps with g(1)>=0 imply phi>=0 and monotonic radial ratios h(xv)/h(x) for x,v>=1. Derivative jumps are harmless: integrate the almost-everywhere derivative over each interval.

The all-scale cyclic reduction is valid with signed angular values. For squared lengths x<=xy<=xz, Q(t)<0 is the angular value at the largest angle; the other two values obey the single sign transition. If the middle value Q(c1)>=0, both later values are nonnegative and their normalized radial multipliers decrease with x. If Q(c1)<0, normalize instead by h(xy)h(xz). The resulting sum is Q(t)h(x)/h(xz)+Q(c1)h(x)/h(xy)+Q(c2). Both ratios increase, their coefficients are negative, and the sum can only decrease. This is a simple alternative to the template's equivalent mixed-sign bracket argument.

For u>=4, c1<=1/2 when t>=-1/4, leaving only the final positive contribution. When t<=-1/4, w²-u²>=3. The final derivative lower gamma=A+min(B0,0)/16 controls the exponential ratios since both arguments exceed16. Thus the exact two tail inequalities suffice; individual coefficient positivity is never needed.

The interval geometry clamp Y,Z>=1 is valid on the entire rectangle: U>=1 and T<=1/2 imply Y=U²>=1 and Z=1+U(U-1)+U(1-2T)>=1. Analytic gradient enclosures are therefore not derivatives of an artificially clipped extension. Cosine clipping to[-1,1] is valid because every such point describes an actual triangle.

The PSD constant block and diagonal correction use C/B² and C/B consistently. Summed rooted forms, paired with Poisson summation, yield P(0)+C plus the pair-cut allowance. The added term373lambda uses the published ordered average count, not half the count and not a pointwise integer bound. Its used half-open shell1.003012 lies strictly inside the source shell1.003013. The outer endpoint is explicitly checked without the allowance. The analytic source theorem is an external dependency; its original SDP is not replayed here.

For squared separation9/8, ball radius sqrt(9/8)/2 has volume pi^(9/2)/Gamma(11/2)*(9/32)^(9/2). The radial function is Schwartz, while the finite angular-hinge kernel has Gaussian decay beyond its knots, so the periodic sums converge. Standard periodic approximation gives the general packing conclusion. None of this proves optimality or worldwide priority.

Actual completion additionally requires compact geometry, Fourier/sign/pair verification and outward rounding on the same frozen files. INPUT_ALIGNMENT.json binds geometry/radial/independent input copies; VERIFIED.json will record the independent compact replay when complete.
