# Audited strengthened near-frame dichotomy: 4/15

Let Λ be a minimum-one lattice in dimension16, with reciprocal covolumeρ<5248 (the existing verified lattice-packing bound suffices). Then

    either ρ<=4096,
    or every16 distinct opposite pairs in a single class modulo2Λ
       satisfy sum_i |v_i|² >=32+4/15.

There may be one representative with squared length below2. This is a proved structural dichotomy, not a universal density improvement by itself.

## Proof

Set G=VᵀV=2I+E, a_i=G_ii−2, S=sum a_i, Q=sum a_i². The same-class midpoint condition implies a_i+a_j>=0 and |E_ij|<=(a_i+a_j)/2. There is at most one negative a_i=−z. In that case S>=14z and Q<=(S−13z)²+15z²<=S²; otherwise Q<=S² directly.

For C=(a1ᵀ+1aᵀ)/2, xᵀEx<=|x|ᵀC|x|, so λ_max(E)<=(S+4sqrt Q)/2<=5S/2. Also tr E²<=8Q+S²/2<=17S²/2 and ||E||<=3S. For S<4/15, G is positive definite. The scalar logarithm estimate gives

    log det(G/2) >= S/2−17S²/[16(1−3S/2)] >=0,

because S<4/15<8/29. Hence detG>=2^16.

The full-rank sublattice L0=V D16/2 lies inΛ and has covolume>=1/128. Its index I is an integer and

    I <= (ρ/128)(1+S/32)^8
       <41*(121/120)^8<44.

Assume S<4/15. An order-p element ofΛ/L0 maps to a nonzero coset of standard D16. For p=11, exact residue histogram enumeration guarantees an unscaled squared-distance witness<=16/11, hence normalized witness<=8/11. For every prime13<=p<=43, the exact parity-support branch certificates guarantee an unscaled witness<=3/2, hence normalized witness<=3/4. The metric distortion from D16/sqrt2 toL0 is at most1+5S/4<4/3. These witnesses would have physical length<1, so none of those primes can divideI.

The remaining torsion cases use exact metric LP dual certificates, which incorporate all midpoint pair inequalities and allow one short representative. Their uniform necessary trace bounds are:

    order5: S>=67/81;
    order7: S>=29/30;
    order9: S>=11/6;
    rank-two3-torsion: S>=186/73.

All exceed4/15. Thus I has no prime factor>=5, no element of order9, and no rank-two3-torsion. Therefore I=2^a or3*2^a. Since I<44, I<=32. Since covolL0>=1/128, ρ=I/covolL0<=4096. This proves the dichotomy.

## Exact computational dependencies

* `verify_metric_families.py five|seven|nine|ternary` reconstructs every rational inequality and checks stored dual certificates without an optimizer, floating-point arithmetic, NumPy or SciPy. There are respectively152,968,4828,4780 cases. Certificate files are `five_metric.json`, `seven_metric.json`, `nine_metric.json`, `ternary_metric.json`.
* `verify_nearframe.py` includes the exact p11 histogram enumeration and earlier conservative torsion checks.
* `verify_prime_support_trees.py` independently regenerates all parity supports and verifies rational branch-tree bounds for p13..43. The geometry and every prime were independently audited. Its companion `prime_support_branch_p*.json` files are necessary proof data.

The large individual metric-LP thresholds MUST NOT be extrapolated as a universal near-frame threshold; determinant and index bounds are separately checked only in the stated range.
