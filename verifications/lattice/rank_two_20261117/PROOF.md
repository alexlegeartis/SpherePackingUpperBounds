# Universal rank-two sign and density-transfer proof

Set h(t)=t^96 exp[-120(t-1)]. Let Q be the selected frozen even polynomial with constant coefficient reduced by2/10000, without normalization. Its Gegenbauer coefficients are nonnegative, q=Q(1)<4e8, and exact Sturm checks give Q<=-eta on[-1/2,1/2], eta=1e-4, and Q'>0 on[1/2,1].

## Scale reduction
A Gauss-reduced basis has Gram w[[1,v],[v,u]], w>=1,u>=1,0<=v<=1/2. Put t=u+1-2v and c0=v/sqrt(u),c1=(1-v)/sqrt(t),c2=(u-v)/sqrt(ut). Main three unoriented lines are e1,e2,e2-e1. Their signed ordered independent contribution is8 h(w)h(wu) times
 Q(c0)+Q(c1)(t/u)^96 e^-120w(t-u)+Q(c2)t^96 e^-120w(t-1).

For w>=1 set A=Q(c1)h(t)/h(u), B=Q(c2)h(t), delta=t-u. If B>=0, the two terms are A p+B q with0<=q<=p<=1, hence <=max(0,A+B). If B<0, then Q(c1)<0 too: c1<=c2, c2>=1/2, Qnegative below1/2 and increasing above. Thus both termsnegative. Consequently Q(c0)<0 and the signed main expression atw1<=-margin imply the same negative sign for allw; quantitative margin is min(eta,margin).

## Omitted-pair tail for1<=u<=9
Every lattice vector outside the6rootvectors either is axial with squarednorm>=4 or is offaxis with squarednorm>=u+1. For offaxis l=+-1 this follows bychecking integerk; |l|>=2 gives norm>=l²(u-v²)>=4u-1>=u+2. Independent pairs between the6rootvectors are exactly the3rootlinepairs. Thus every omitted independent pair has normalized energy E>=E0+1, E0=1+u<=10.

One of the two vectors is offaxis, so their squarednorms can be ordered r>=1,s>=u. Put Delta=E-E0. Concavity/log-AMGM gives
 h(wr)h(ws)/(h(w)h(wu)) <= (1+Delta/2)^192 exp(-120wDelta).
For Delta>=1 this is <=R exp[-56(Delta-1)], R=(3/2)^192 exp(-120).

Disjoint disks of radius1/2 around lattice points show N(norm²<=T)<=9T for T>=1; ordered-pair counting by energy gives N_pair(E)<=81E². Stieltjes integration therefore bounds omitted absolute positive contribution/base by
 81 q R[(E0+1)²+2(E0+1)/56+2/56²]
 <=81 q R[121+22/56+2/56²]<2e-6.
Thus main signed expression<=-1e-6 contributes<=-8e-6 and dominates.

For3<=u<=9, c0,c1<=1/2, so the main expression is at most -eta+q h(3), certainly below-1e-6. (h decreases for t>=1.) Only1<=u<=3 requires a compact interval gate; the included interval gates cover1<=u<=4.

## Anisotropic u>=9
A positive-Q pair requires |cos|>1/2. If both vectors offaxis, normalized energy>=2u. If one axial, its partner's norm>=4/3(u-1/4), so energy>=1+u+(u-1)/3. Therefore positive-pair energy starts at Emin=E0+Delta, Delta=(u-1)/3. log h'<=-24 for t>=1 gives weight/base<=exp[-24w(E-E0)]. Counting and Stieltjes integration bound positive sum/base by
 81q exp[-24wDelta][Emin²+Emin/12+1/288].
At w>=1,u>=9 this decreases inw,u, and atw1,u9 it is<1e-15. Atleast8 signed ordered shortest-axis/secondbasis pairs contribute<=-8eta. Allothersnonpositive can be dropped. This proves negativity in the anisotropic region.


# Dimension-adapted kernels and strengthened origin transfer

The same universal rank-two theorem applies to each frozen polynomial after subtracting2/10000 from its Gegenbauer constant coefficient. Q45 uses Gegenbauer parameter43/2; Q42 uses parameter20. Each polynomial's positive-definite expansion in a target dimension is checked exactly, and it is then normalized by its target-dimension constant coefficient. The raw geometric sign proof uses q<4e8, eta=1e-4 and unchanged h(t)=t^96 exp[-120(t-1)].

Frozen Q45 SHA256: defe0aa8ab8d91f7bd14a983c341bf749c3695d01af509431c5ae77e0b6f1712.
Frozen Q42 SHA256: 3d9ee200d60ab7784a169bf73059a4e7e0688423060f28987121d1828d2162cf.
Both were recovered after a computer rollback by deterministic reoptimization, with exact byte-hash matches to independent audits. Restored inputs are not merely numerical approximations to the originals.

## General collinear domination
Let old scalar kernel parameters be c=C/B², u0=c*qold, v0=c*R, where R is the actual saved origin ratio timesB. The prior per-line correction is
 4u0 sum h_k² -4v0 sum h_k +4cB sum h_k h_(2k).
For new normalized polynomial Q with qnew=Q(1), use a=u/qnew, cross feature -v*h, and origin s=qnew*v²/u. The new complete line correction is
 4u (sum h_k)² -4v sum h_k.
Set rho=1/100 and Gamma=(1-rho)(1-rho²). The established geometric decay gives sum_(k<l) h_k h_l <= h1*h2/Gamma.

If u>=u0, v-v0=u-u0, and2u<=cB*Gamma, the new full correction is dominated by the old extracted correction: their difference equals
 4(u-u0) sum(h_k²-h_k) +8u sum_(k<l)h_k h_l -4cB sum h_k h_(2k),
and both combined terms are nonpositive. Here0<h_k<=1. This is a global line-sum inequality, exactly what the PSD argument needs.

For ordinary invariant transfers, chooseu=u0,v=v0. The necessary sufficient check is2qold/B<=Gamma, not the earlier convenient but unnecessarily strong3/4 threshold.

For40D, old origin ratio is1 andu0<cB*Gamma/2. Chooseu=cB*Gamma/2 andv=v0+(u-u0). This decreases the origin beyond the invariant transfer, while preserving domination of the same radial/profile witness. check_rank2_transfer40_optimized.py verifies all parameters and rational safe rounding.


## From rank-two signs to a lattice density bound
For any lattice Λ of minimum at least1, partition every linearly independent ordered pair by its real two-dimensional span. Each saturated intersection with that span is a rank-two lattice of minimum at least1. Absolute Gaussian convergence justifies regrouping. The universal sign theorem makes the sum of all independent-pair kernel terms nonpositive.

A nonnegative Gegenbauer expansion makes K(x,y)=a h(|x|²)h(|y|²)Q(cos angle) positive semidefinite. Choose its origin cross term −v h and origin scalar s=q v²/u, where a=u/q and q=Q(1) after normalizing the target-dimensional constant coefficient to1. The constant-harmonic2×2 origin block is PSD. Apply PSD to growing finite lattice sets and pass to the absolutely convergent sum. After dropping the nonpositive independent contribution, the origin plus full primitive-line correction is nonnegative. The preceding domination lemma compares this correction with the exact old radial majorant.

Poisson summation for the frozen Schwartz radial Fourier pair gives sum f(Λ)≥fhat(0)/detΛ. Combining the PSD inequality and verified mod2 class-profile budgets bounds reciprocal covolume, hence packing density. Each predecessor radial and profile witness is rechecked; no predecessor optimality or density estimate is assumed as a theorem premise. Its previously certified safe decimal is used conservatively for the exact objective-ratio conversion, including odd dimensions.

## Direct36D Fourier/profile refinement
Dimension36 uses the Q37 rank-two sign theorem restricted exactly to36D, but a newly verified quadratic nonnegative Fourier multiplier and radial/profile certificate, not invariant origin transfer. Its complete original four angular-cut inequalities are replayed, then combined with the new nonnegative cut weights. The class budget includes2 sum ν_i constant_i separately from2(2^36−1)κ. The separate exact Laguerre inversion checks the new Fourier transform. See multiplier36_independent/verify.py and check_fourier_transform.py.

## Verification status and scope
All five scalar angular gates and required predecessor radial/profile dependencies passed the independent audits preserved in the source backups. The release runner repeats the complete proof chain. These scalar kernels are distinct from subsequent harmonic-dependent matrix candidates; no matrix candidate is included. Comparisons are to the explicitly pinned preceding repository certificates, not assertions of exhaustive worldwide priority.

##46D and48D Fourier refinements
The46D quadratic multiplier is PSD and nonnegative on the Fourier domain; a separate exact Laguerre inversion checks the physical transform.48D adds a polynomial of the form b4(s−1)^T M b4(s−1)+s b3(s−1)^T N b3(s−1), with exact positive-definite5×5 and4×4 matrices. The quadratic multiplier constant is adjusted so the total Fourier value at zero is unchanged. Both replays check global Fourier positivity, the new physical tail, every finite radial bin, origin PSD and all class-profile duals. These are direct certificates, not merely the old radial origin transfer. The universal scalar rank-two theorem and full-line budget remain unchanged. See the included independent audit notes and exact transform checks.
