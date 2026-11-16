# Global short-class angular coupling with exact missing-target penalties

This is a rigorous structural inequality. Numerical polynomial searches using it are separate and require continuous sign certificates before they establish packing bounds.

Let Λ⊂R^n have minimum length1. Fix1≤U<2 and let A be the set of nonzero classes modulo2Λ whose unique short opposite pair has squared norm in [1,U). Write N=|A| and choose one vectorv_a per pair. Each nonzero class contains at most one opposite pair with squared norm<2.

Let F(t)=Σ_{k even} b_k P_k^(n)(t), with b_k≥0, where P_k is normalized Gegenbauer, P_k(1)=1. Put b0=b_0 and H=F(1). Assume

    F(t)≤0 for0≤t≤U−1.

Partition the target short classes by their squared norms into disjoint half-open bins I_j=[l_j,u_j)⊂[1,2). For each bin choose W_j≥0 bounding the positive contribution from ordered pairs(a,b)∈A²,a≠b,a+b=c, uniformly for each targetclassc in that bin. Then

    b0 N² ≤ H N + Σ_j W_j M_j,                 (1)

where M_j is the number of target short classes in binj. Classes without a short pair contribute no positive penalty.

## Proof of no penalty for nondependent/empty targets

For distincta,b set t_a=|v_a|²,t_b=|v_b|² and s=min|v_a±v_b|²=t_a+t_b−2|<v_a,v_b>|. If s≥2, then

    |<v_a/|v_a|,v_b/|v_b|>|≤(U−1)/U≤U−1,

so F contributes nonpositively. If s<2, the vector attaining s is the unique short pair in targetclassc=a+b. In particular a targetclass with no short pair cannot contribute positively. Every positive contribution therefore comes from an actual equalityw=u+v after uniquely orienting sourcevectorsu,v and fixing an orientationw of the targetpair. Merely knowing a+b=c would not suffice; the positivity/sign argument is essential.

The set of ordered sourceclass pairs contributing to a fixed orientedw injects into decompositionpoints u with u,v=w−u bothsourcevectors. Conversely each suchu identifies its orderedclasses. Swappingtheclasses gives the other decompositionpointv; there is NO extra factor2 beyond countingbothpoints. Summing the Gegenbauer positive-definite inequalities on one representative per sourcepair yields b0 N²≤HN+Σ_{a≠b}F(cosangle), proving(1) after bounding its positive terms.

## Angular interval for a dependent target

When |w|²∈[l_j,u_j], |u|²,|v|²∈[1,U], a safe interval for the absolute normalized inner product is

    J_j=[max(0,(2−u_j)/(2U)), min(1−1/(2U),(2U−l_j)/2)].

Choose B_j≥max(0,sup_{t∈J_j}F(t)). If B_j=0 takeW_j=0. Otherwise multiplyB_j by any of the exact multiplicitycaps below.

## Basic combinatorial cap

For each a there is at mostone b=a+c, so at mostN orderedpairs use a fixedtarget. This gives the weaker inequality b0N²≤HN+NΣB_j M_j, hence b0N≤H+ΣB_j M_j whenN>0.

## Antipodal ball cap

For fixedw the decompositionpoints lie in the ball centeredw/2 with radius squaredR²=U−|w|²/4≤U−l_j/4. They are centrallysymmetric underu↔w−u and mutuallyseparated by at least1. Their centeredantipodalpairs ±y have |y|≥1/2. For two distinctpairs, |y_i±y_k|≥1, so their normalized directions have absoluteinnerproduct at most

    c_j=1−1/(2R_j²),  R_j²=U−l_j/4,

providedR_j²<1. Monotonicity of(a+b−1)/(2sqrt(ab)) in both squaredradii follows froma−b+1>0 whenb<1. If0≤c_j andn c_j²<1, the projectiverelativebound gives

    C_j=2 floor[n(1−c_j²)/(1−n c_j²)]

as an upper bound on the NUMBEROFPOINTS, hence orderedpairs. IfR_j²<1/2 there can be at mostoneantipodalpair, givingC_j=2. ThusW_j=C_j B_j is valid. No sampled spherical-code optimization is needed for these caps.

## Stronger cap using the existing conditional mod2 class budget

For each decomposition defineq=2u−w=u−v. Thenq lies in the SAME mod2classasw, q≠±w (otherwiseuorv=0), and

    |q|²=2|u|²+2|v|²−|w|²≤4U−l_j.

Oppositeq pairs correspond to two ordereddecompositions. Consequently one can instead take

    C_j=2 A_n(sqrt(u_j),sqrt(4U−l_j)),

where A_n is the standard conditional classcap (known shortpair excluded). Use disjoint half-open squared-norm bins [E_j,E_{j+1}) for class accounting. The short part is [1,2), the long part is [2,E_end), and the nonpositive tail covers E_end itself. Actual short and long squared norms are strictly below their respective bin upper endpoints. Thus the conditional zero-cap test on endpoint bounds a+S<=4 is valid: the actual squared norms have sum strictly less than 4. At actual equality 4, distinct orthogonal pairs can occur, so the test would not be valid for closed-bin counts. Angular and radial polynomial inequalities are checked on the closed containing intervals, which is conservative. Likewise a source-prefix window denotes [1,E_prefix); the recorded closed source_interval is an envelope, not a second assignment of boundary classes. Its projective projection formula applies whenever its denominator ispositive; if it fails, use the basic cap or another certified bound rather than asserting a finite relativebound. This is the direct bridge between global class addition and the existing within-class cap.

## Linear cuts for the profile LP

For any fixedN0≥0, (N−N0)²≥0 and(1) imply

    (2b0N0−H)N−Σ_jW_jM_j≤b0N0².            (2)

Takingb0=1,N0=H gives N≤H+Σ_j(W_j/H)M_j. More generally fixN0 near the currentprofile. The coefficients ofF and perbinupperboundsB_j can be jointly optimized with the radialcertificate: b_k≥0, F≤0on[0,U−1], andB_j≥F onJ_j are linear constraints. The capsC_j are fixedrational/integerdata. Thus(2) is linear in the newdualvariables; their scale absorbs any multiplier.

For a sampled search, all angularinequalities remain exploratory. A final certificate mustverify F≤0on[0,U−1], eachB_jmajorant onthefullclosedintervalJ_j, Gegenbauercoefficientnonnegativity, and all relativebounddenominators exactly. Countnormalizations mustuseone sourceoppositepair perclass, orderedAAoffdiagonalpairs, and diagonaltermHN.

## Sharper rational cap by projection perpendicular to the target

An independent audit gives a slightly stronger version of the preceding projection cap without separate endpointslack. Write y=u−w/2 and t=|w|². Then

    |y·w| = ||u|²−|w−u|²|/2 ≤(U−1)/2,
    |y|²=(|u|²+|w−u|²)/2−t/4≥1−t/4.

Thus the squared directional correlation with w is at most

    ε_j=(U−1)²/[l_j(4−l_j)],

because t(4−t) is increasing on[1,2]. If s_j=1−1/[2(U−l_j/4)], the normalized projections of the antipodal directions onto w^perp have absoluteinnerproduct at most

    c'_j=(s_j+ε_j)/(1−ε_j).

When0≤s_j, ε_j<1 and(n−1)(c'_j)²<1, use

    C_j=2 floor[(n−1)(1−(c'_j)²)/(1−(n−1)(c'_j)²)].

The projection is injective on distinctantipodalpairs because the resulting absoluteinnerproduct bound isstrictlybelow1. All quantities are rational for rationalbinendpoints andU. This incorporates the narrowwidth of the sourcebin and is generallystronger than the unprojectedn-dimensional relativebound.

## Dimension-transfer family of explicit projective caps

The Welch and quartic formulas are the first two members of a general family. For integer r≥0, set

    Q_r(x)=x^(2r)(x²−c²).

It is nonpositive on[−c,c]. In dimensiond≥3, the coefficient of normalizedGegenbauer P_(2j) in x^(2r+2), divided by its coefficient in x^(2r), is

    (r+1)(2r+1)/[(r+1−j)(d+2r+2j)],  0≤j≤r.

The denominator differs from its j=0 value by j(2−d)−2j²≤0. Therefore the ratio is minimizedatj=0 and every Gegenbauercoefficient ofQ_r is nonnegative whenever c²<(2r+1)/(d+2r). The new top-degree coefficient ispositive. Its constantcoefficient is

    μ_(2r) [(2r+1)/(d+2r)−c²],
    μ_(2r)=(2r−1)!!/[d(d+2)...(d+2r−2)],   μ_0=1.

The projectivepaircount is consequently bounded by

    K_r=(1−c²)/{ μ_(2r)[(2r+1)/(d+2r)−c²] }.

Use C=2 floor(min_r K_r) for the ordereddecompositioncount. Forr=2 this is the sexticformula

    C=2 floor[d(d+2)(d+4)(1−c²)/(3(5−(d+4)c²))].

At sourceU=97/96 and worsttargetlowerendpoint1, the projectedc is exactly691273/2018231. The sexticformula covers everytargetbin throughambientdimension39; higher r extends the coverage further. The packaged complete verifier reconstructs the required quadratic, quartic and sextic Gegenbauer coefficients exactly. These are standard Delsartepolynomial constructions, not claims that the formula family itself is novel. Their use here is to give portable exactcaps in the globalclasscoupling.
