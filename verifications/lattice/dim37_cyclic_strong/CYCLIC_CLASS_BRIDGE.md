# Audited bridge: cyclic PSD correction plus lattice class budgets

This is a theorem-level derivation, not a completed numerical certificate.
Let Lambda be a Euclidean lattice of minimum squared norm1, reciprocal
covolume rho, and m=2^n-1. Let f be a real even Schwartz function with
Fourier(f)>=0 and Fourier(f)(0)>0.

Assume K is a real symmetric origin-augmented PSD kernel with K(0,0)=C>=0,
absolute convergence of all following lattice sums, and
 K(x,y)+K(-x,y-x)+K(-y,x-y)<=0
whenever x,y,x-y are nonzero and each has squared norm>=1. This includes
collinear triangles. Put D(v)=K(v,v)+2K(0,v).

For the ordered set of pairs x,y in Lambda\{0}, x!=y, cyclic changes of anchor
are bijections. Symmetry of K identifies the third ordering correctly. Summing
the cyclic inequality thus gives O=sum_{x!=y}K(x,y)<=0. Finite PSD sums,
followed by absolute-convergence limits, give C+sum_{v!=0}D(v)+O>=0.
Consequently sum_{v!=0}D(v)>=-C. No averaging across centers, finite-window
boundary assumption, or individual off-diagonal sign condition is needed.

Suppose G is an even radial majorant f(v)+D(v)<=G(|v|^2) on |v|^2>=1, and
valid lattice class-profile inequalities imply sum_{v!=0}G(|v|^2)<=2m*kappa.
Poisson summation and nonnegative Fourier transform yield
 rho*Fourier(f)(0) <= sum_v f(v)
 <= f(0)+C+sum_{v!=0}G(|v|^2)
 <= f(0)+C+2m*kappa.
Each nonzero vector is counted once; the factor2 in 2m accounts only for the
opposite pairs used in class budgets. There is NO additional factor2 on C.

## Published profile normalization
For Q=sum c_l G_l^(n), c_0=1,c_l>=0, B=Q(1), and h(v), the published kernel is
 K(x,y)=C/B^2*h(|x|^2)*h(|y|^2)*Q(angle),
 K(0,x)=-C/B*h(|x|^2), K(0,0)=C.
Thus K(x,x)=C/B*h^2 and D=C/B*(h^2-2h). The constant-harmonic feature block
[[C,-C/B],[-C/B,C/B^2]] is PSD of rank one. Other harmonics are PSD by the
addition formula. For h(v)=v^96 exp[-120(v-1)], h decreases from1 on v>=1,
so D<=0 there and has Gaussian decay. If f<=0 past a cutoff, f+D<=0 there
also. The published cyclic triangle gate, including all dilations and the
infinite tail, must be replayed or explicitly retained as a verified dependency
for the IDENTICAL dimension, Q and h. It does not transfer to a different
angular dimension or changed radial profile without a new proof.

If optimizing variable C, the spatial column is (h^2-2h)/B and its objective
coefficient is1. If optimizing A=C/B, the spatial column is h^2-2h and the
objective coefficient is B. Coefficients must remain nonnegative.

The G=0 choice contains the published general-packing certificate when its
radial polynomial degree and scale are retained. A smaller degree or a fixed
rational Gaussian scale can exclude this control; after minimum-distance
normalization the original scale may be 2*pi*s. Therefore worse numerical
results under those restrictions do not refute the bridge's usefulness.
