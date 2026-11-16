# Dimension-uniform near-frame determinant lemma

## Attribution and novelty correction (2026-11-16)

The determinant step is already contained in Rémi Géraud-Stewart and David Naccache, *A Focus on Hermite Constant 11: the certificates and macadamia 0052* (8 November 2026), Appendix D, Lemma 32, pp. 44–45: https://eprint.iacr.asia/2026/2826.pdf . For n >= 7, the contact inequalities E_ii+E_jj >= 2|E_ij| and s=tr(E) in [0,1] imply det(I+E) >= 1+s-(n-1)s^2/4 when the right side is positive. Substituting E=G/2-I and s=S/2 gives the determinant range used here. Its proof already uses positive/negative-spike diagonals, fixed off-diagonal ratios, and log-determinant concavity. The derivations/checkers below are independent reconstructions, not a novelty claim for this lemma or proof. Any potential novelty must lie elsewhere and needs a separate literature comparison.

## Lemma

Let n≥9, and let G be a real symmetric n×n matrix with

    Gii=1+δi/2,
    |Gij|≤(δi+δj)/4  (i≠j),
    δi+δj≥0,  Σδi=S,  0≤S≤8/(n−1).

Then G is positive definite and detG≥1. Individual δi may be negative.

## Convex reduction

When S=0, pairwise nonnegative sums imply allδi=0 (n≥3), so G=I.

For S>0, the δ-polytope is the convex hull of

    A_i=S e_i,
    B_i=S/(n−2)(1−2e_i).

If allδi≥0, use the A_i combination. Otherwise let δi=−a be its unique negative entry. All otherδj≥a and (n−2)a≤S. Use weight λ=(n−2)a/S on B_i, followed by weights (δj−a)/S on A_j for j≠i. These weights are nonnegative and sum1.

For each pair choose θij∈[−1,1] such that Gij=θij(δi+δj)/4; if the denominator vanishes then Gij=0 and chooseθij=0. Keep these ratios fixed during the convex decomposition. Then G(δ) is affine inδ. It suffices to prove positive definiteness and determinant≥1 for every extreme matrix, because logdet is concave on positive definite matrices.

## Positive-spike extreme

For δ=A_i, G has an identity (n−1)×(n−1) principal block and a star row b with |bj|≤S/4. Its Schur complement is

    1+S/2−||b||² ≥1+S/2−(n−1)S²/16 ≥1.

Hence it is positive definite with determinant≥1.

## Negative-spike extreme

Put t=S/[2(n−2)], so 0≤t≤b=4/[(n−1)(n−2)]. For δ=B_i, the distinguished coordinate decouples and has diagonal1−t. The remaining m=n−1 block H has diagonal1+t and off-diagonal absolute values≤t.

Write H=(1+t)(I+xA), x=t/(1+t), where A has diagonal0 and |Aij|≤1. Let J be the matrix with diagonal0 and all off-diagonal entries1. Then |tr(A^k)|≤tr(J^k) by the closed-walk expansion. Since x(m−1)<1,

    log det(I+xA)
      =Σ_{k≥2} (−1)^(k+1) x^k tr(A^k)/k
      ≥−Σ_{k≥2} x^k tr(J^k)/k
      =log det(I−xJ).

Positive definiteness follows independently by Gershgorin. Consequently

    det H ≥(1−(n−3)t)(1+2t)^(n−2),
    det G ≥D_n(t):=(1−t)(1−(n−3)t)(1+2t)^(n−2).

All factors are positive: (n−3)b<1 for n≥9. We show D_n is increasing on[0,b], hence D_n(t)≥D_n(0)=1.

The numerator of d(log D_n)/dt, over the positive denominator (1−t)(1−(n−3)t)(1+2t), is

    N_n(t)=(n−2)−2(n²−4n+5)t+2n(n−3)t².

Its derivative is increasing in t, and at t=b it equals

    −2(n^4−7n^3+11n²+n+10)/[(n−2)(n−1)]<0.

Thus N_n decreases on[0,b]. Its endpoint value is

    N_n(b)=(n^5−16n^4+81n^3−158n²+116n−88)/[(n−2)²(n−1)²]>0.

Both polynomial signs hold uniformly for n≥9. Writing n=k+9, k≥0, their numerators become respectively

    k^4+29k³+308k²+1414k+2368,
    k^5+29k^4+315k³+1543k²+3104k+1280.

Every coefficient is positive. This proves the lemma.

## Lattice interpretation

For n distinct opposite pairs v_i in one mod-2 lattice class, writeδi=|v_i|²−2. Then T=[v_1 … v_n]/sqrt2 has TᵀT of the stated form. If S≤8/(n−1), T is invertible and detT≥1. Moreover T(D_n/sqrt2)⊂Λ and its covolume is

    2^(1−n/2) detT ≥2^(1−n/2).

This determinant lemma is not by itself a packing-density improvement; it supplies the dimension-uniform geometric part of a torsion/index exclusion.
