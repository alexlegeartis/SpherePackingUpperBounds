# Independent antipodal allocation under the basic count gauge

This is a local-energy construction and a numerical dual-pricing diagnostic. It is not a new packing bound. The frozen geometry used below comes from a finite numerical master; a useful reduced cost gives a direction to test, not an achievable density improvement.

## The local-energy identity

Let (C) be the nonzero neighbors of a selected packing center, restricted to the active annulus (1\leq |x|\leq R). Their mutual distances are at least one. Put

\[
\Phi(z)=(1-2|z|^2)_+.
\]

The elementary Blichfeldt capacity inequality gives
\(\sum_{x\in C}\Phi(x+y)\leq1\) for every (y). Therefore, for every nonnegative radial function (u) on the annulus,

\[
E_u(C)=\sum_{y\in C}u(|y|)
       \left(1-\sum_{x\in C}\Phi(x+y)\right)\geq0.
\]

The origin constant is zero, the origin-to-neighbor term is (h_u(y)=u(|y|)/2), and the symmetric neighbor kernel is

\[
K_u(x,y)=-\frac{u(|x|)+u(|y|)}2\Phi(x+y).
\]

The kernel is defined by this formula only when both nonzero arguments lie in the active annulus, and it is zero otherwise. Since \(|y|\geq1\), the self value \(\Phi(2y)\) is zero. Thus adding a multiple \(\gamma E_u\), with \(\gamma\geq0\), raises the complete unlifted pair expression by exactly \(\gamma u(r)\).

Retain the basic count capacity \(\psi(x)=\mathbf1_{[1,R_1]}(|x|)/24\), and write \(d(r)=1-\psi(r)>0\). If the original complete unlifted pair expression is \(q(r)\), the new retained pair inequality is

\[
q(r)+\gamma u(r)\leq0.
\]

Its fresh count-gauge weight is

\[
w_\gamma(r)=-\frac{q(r)+\gamma u(r)}{d(r)}
            =w_0(r)-\gamma\frac{u(r)}{d(r)}\geq0.
\]

Consequently the actual lifted pair remains identically zero. For a neighbor pair with radii \(r_i,r_j\), the change in the actual kernel per unit \(\gamma\) is

\[
\delta K_{ij}=\frac12\left(
 \frac{u_i\psi_j}{d_i}+\frac{u_j\psi_i}{d_j}
 -(u_i+u_j)\Phi(x_i+x_j)\right).
\]

There is no origin cost. The positive terms in this display are the count budget relinquished to pay for the independent allocation. They must be included, as must the new unlifted pair row (u).

## Full triangles, exterior triangles, and lattice rows

At the present support \(R=7/4\), a full triangle satisfies

\[
|x_i+x_j|^2=2r_i^2+2r_j^2-r_k^2
           \geq4-R^2=15/16>1/2.
\]

Hence the antipodal correction vanishes on every full triangle, including all boundary closures. The full cyclic change is solely
\(\sum_i u_i\sum_{j\ne i}\psi_j/(2d_i)\).

On an exterior branch parameterized by

\[
\cos\theta=-1+z\frac{(r+s)^2-R^2}{2rs},\qquad0\leq z\leq1,
\]

the exact hinge is

\[
\Phi(x+y)=\left[1-2(r-s)^2-2z\big((r+s)^2-R^2\big)\right]_+.
\]

This formula also covers the antipodal face (z=0), without a division by an angular or radial gap.

A radially truncated lattice neighborhood is centrally symmetric. For every selected (y), the target (-y) contributes one to the capacity, and every other target contributes zero because its distance from (-y) is at least one. The independent energy (E_u) is therefore exactly zero on these neighborhoods. The only lattice-row change is the lost old count-gauge energy,

\[
-\left(1-\sum_{x\in C}\psi(x)\right)
  \sum_{y\in C}\frac{u(|y|)}{d(|y|)}.
\]

The numerical script evaluates this using the existing shell counts and exact row-region labels. It does not reconstruct side conventions from rounded point coordinates.

## Frozen-dual result

The script `price_independent_antipodal.py` prices the ten old quartic Bernstein functions and three old count indicators against `nested_conditional_jointcount_fullgauge_dual_geometry.npz`, SHA-256 `9608d101717d76bb61a4b2ece9bb80c61598c69886cbc354dab9cd0266c9ac27`.

Exactly one of these thirteen profiles has a negative reduced cost: the fourth-degree endpoint Bernstein function on the first broad interval,

\[
u(r)=(2r-2)^4\quad(1\leq r\leq3/2),
\]

with the old first-interval boundary convention, and zero on the second interval. It is the profile labeled `inner_B4`.

| Component | Raw reduced cost |
|---|---:|
| Retained unlifted pair row | +0.014352418782534 |
| Full-triangle count-budget loss | +5.330794792597345 |
| Exterior count-budget loss | +5.267254211499037 |
| Exterior antipodal energy | −11.435136941123849 |
| Lattice count-budget loss on this dual | 0 |
| Total | **−0.822735518244933** |

The raw objective is eight times the density objective, so the corresponding first derivative is approximately (-0.102841939780617) per unit profile weight. This is a local numerical price; the admissible coefficient is limited by the fresh pair inequalities and all other primal inequalities.

All twelve other tested profiles price positively. The full result is saved in `independent_antipodal_allocation_prices.json`. The computation covers 1,942 pair rows, 9,606 full rows, 16,402 exterior rows, and 280 lattice rows. Reconstructing the independently exported physical origin gradient from the pair, count-budget, and ordinary lattice terms agrees to (1.42\times10^{-14}). This verifies the numerical sign and factor-of-two conventions. It does not replace exact interval certification of any future candidate.
