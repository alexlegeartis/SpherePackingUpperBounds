# Compact harmonic-matrix lattice correction in48D

Producer-complete; independent audit pending. Proposed safe upper3.499004606e−7. The earlier coarse sampled matrices failed adversarial contact tests and are not used. frozen.json is the repaired final input.

Let R=19/13, h(t)=max(1−13(t−1)/6,0)^12 for t>=1, and H(t)=(h(t),h(t)^2). Let c_k be the nonnegative normalized even Gegenbauer coefficients in frozen.json, c_0=1, q=sum c_k. All symmetric2x2 matrices B_k are positive definite by exact rational principal-minor checks. The exact origin block [[s,z^T],[z,B_0]] is positive semidefinite by its rational Schur complement.

Define K(x,y)=q^−2 sum_k c_k G_(2k)(cos(x,y)) H(|x|²)^T B_k H(|y|²), K(0,x)=q^−1 z^T H(|x|²), K(0,0)=s. Matrix-valued Gegenbauer positivity gives a positive semidefinite kernel. All nonzero radial features have compact support squared length<R<2.

## Complete lattice-plane coverage
For a lattice of minimum length1, choose a Gauss-reduced basis in each rank-two plane. Write its squared lengths a<=b and dot product d with0<=d<=a/2; put c=a+b−2d>=b. Every vector other than ±e1,±e2,±(e2−e1) has squared length>=2a>=2, and is therefore unsupported. The proof is the standard reduced-coordinate bound: |l|>=2 gives (b/a−(d/a)^2)l²>=3; |l|=1 outside the two designated roots gives squared length>=a+b>=2; multiples of e1 have squared length>=4a.

If fewer than two root pairs are supported there is no independent contribution. If all three are supported,1<=a<=b<=c<=R. Parameterize the whole domain by x,y,z in[0,1]: a=1+(6/13)x, b=a+(6/13)(1−x)y, c=b+(6/13)(1−x)(1−y)z. Then h(a)=(1−x)^12, h(b)=h(a)(1−y)^12, h(c)=h(b)(1−z)^12. Dividing the independent signed-pair sum by the positive8h(a)h(b)/q² yields exactly parts(...,mode=main) in gate.py. It has complete interval coverage with122609 nodes and61305 leaves.

If exactly two root pairs are supported,1<=a<=b<=R and c>=R. Gauss reduction gives0<=d<=(a+b−R)/2. The full complementary domain is parameterized by the same a,b and d=(a+b−R)z/2, z in[0,1]. The normalized independent-pair contribution is exactly mode=cap. Its complete interval cover has289463 nodes and144732 leaves.

Evenness removes all square roots from angular evaluations. Boundary points with h(a)h(b)=0 are limiting normalized formulas; their actual independent contribution is zero. No proof step divides by zero. Both gates use exact rational endpoints,256-bit Arb bounds, centered polynomial evaluation, and complete first/second-derivative Taylor enclosures. They reject a positive center immediately. Full preorder subdivision traces and input hashes are saved. All independent ordered pairs partition by their lattice plane, so the nonpositive plane sums imply the global independent-pair sum is nonpositive.

## Collinear sum and scalar closure
Compact support<R<2 leaves at most±one supported vector on each lattice line. Put W=q^−1 sum_k c_k B_k. The exact pointwise correction is G(t)=2q^−1(H(t)^T W H(t)+z^T H(t)); no multiple-line tail estimate is needed. Kernel PSD and nonpositive independent sum imply sum_(v!=0)G(|v|²)>=−s.

The repaired Fourier polynomial is nonnegative globally by rational positive factors and has value1 at0. Exact Laguerre inversion constructs f. The radial verifier proves f+G<=height_i on633 annuli, with the exact cutoff x=418/13 included as a partition boundary. Monomial antipodal-code bounds in each mod2 class, counting all2^48 classes conservatively, give total budget B=2^49 sum cap_i height_i. For annular upper squared radius<2 the cap is1: two distinct unoriented vectors in a class would make both half-sum and half-difference nonzero lattice vectors, contradicting parallelogram. All potentially-real roots of the physical polynomial lie below180 and P(180)<0, while G vanishes there, establishing the full infinite tail.

Thus density <=(22/4)^24/24!*(f(0)+s+B). Every scalar and origin coefficient is exact rational, and the final rounding is checked strictly. Additional shape/symmetry assertions were added to the verifiers after the producer run; they do not change the evaluated functions or coverage. Independent replay should include them.
