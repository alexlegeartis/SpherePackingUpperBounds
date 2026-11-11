# Methods, proof dependencies, and next experiments

This document organizes the repository's packing-density methods and a finite research agenda. It does not replace [RESULTS.json](RESULTS.json), the dimension-specific proofs, or `python run_all.py`. Proposed experiments below are unverified; a numerical objective is never a density bound. No optimality, external peer-review, or worldwide-priority claim is made.

## Current 17D–20D certificate chain

| Dimension | Earlier solid-harmonic cyclic bound | Refined strict bound | Refinement and proof location |
|---|---|---|---|
| 17 | 0.018335552 | 0.018204096 | Wider support, explicit mixed-Gaussian tail repair, then exact origin saturation: [17D package](verifications/dim17/refinement_repaired/) |
| 18 | 0.013284728 | 0.013224473 | Exact origin saturation: [common package](verifications/cyclic_origin_slack/), using [18D cyclic source](verifications/dim18/refinement_cyclic/) |
| 19 | 0.009557491 | 0.009521000 | Exact origin saturation: [common package](verifications/cyclic_origin_slack/), using [19D cyclic source](verifications/dim19/refinement_cyclic/) |
| 20 | 0.006977702 | 0.006936042 | Exact origin saturation: [common package](verifications/cyclic_origin_slack/), using [20D cyclic source](verifications/dim20/refinement_cyclic/) |

The earlier certificates remain valid and preserved. In particular, the earlier 17D witness alone admits the saturated bound 0.018296869; the stronger registered refinement uses a different, wider-support witness and its explicit radial repair. These are distinct certificate chains, not interchangeable inputs.

## One proof skeleton, several finite function families

The common construction combines a Fourier-positive radial function f with a positive-semidefinite kernel K. With H0 = fhat(0)>0, local pair inequalities and nonpositive global triangle sums give center density at most

    (f(0) + K(0,0) + count_cost) / H0.

The packing-density bound is this expression times the volume of the ball whose radius is half the minimum center separation. Periodic approximation extends the conclusion to unrestricted packings, rather than only lattices.

The framework is [Cohn, de Laat and Salmon, Three-point bounds for sphere packing](https://arxiv.org/html/2206.15373v1). Its printed Theorem 1.1 requires individual kernel nonpositivity. Section 2.2, equation (2.4), and Theorem 2.4 justify replacing that by the weaker cyclic-sum condition. The repository also gives the finite periodic summation explicitly, including the isolated origin and support boundaries. The cyclic mechanism is existing mathematics; the computational witnesses are not new density theorems.

### 1. Fourier-positive two-point baseline

Setting K=0 recovers the Cohn–Elkies construction. A polynomial times a Gaussian has an exact Laguerre/Fourier transform, and two global ray-sign checks imply the bound. Low-degree Gaussian and sparse reciprocal-Gaussian experiments are useful controls but did not produce competitive bounds in the tested cases. Dense numerical sign screening rejected by a subsequent exact root test is a failed candidate, regardless of solver status.

### 2. Annular pointwise-negative kernels

The [annular search package](research/annular_three_point/) uses piecewise-linear radial hats in squared radius, an isolated degree-zero origin coordinate, and normalized Gegenbauer blocks. Individual angular signs together with bilinear interpolation cover the complete pair domain. This is the earlier 18D, 19D, 20D and 21D family. The radial coupled term is quadratic on each cell. These certificates do not use external near-neighbor count bounds.

Relevant proofs are in [18D](verifications/dim18/), [19D](verifications/dim19/), [20D](verifications/dim20/), and [21D](verifications/dim21/). Historical proof files describe their own frozen witness, not necessarily the newest registered bound in their dimension.

### 3. Solid-harmonic cyclic kernels

The [solid-cyclic search package](research/solid_cyclic/) multiplies angular degree ell by (sqrt(uv)/s)^ell, where u and v are squared radii and s is the exclusion squared. Positive radial feature factors preserve PSD. The resulting angular terms are rational polynomials in squared triangle sides, enabling exact Bernstein certification of cyclic sums.

There are separate proof obligations for compact triangles and triangles with one edge outside the support. With two long edges all vertex terms vanish. The exact 17D witness has a rational test triangle with a positive individual vertex kernel and a negative cyclic sum. This proves that the witness genuinely uses cyclic relaxation; it neither proves optimality nor isolates the fraction of its improvement attributable to the relaxation.

### 4. Count-enhanced low-dimensional kernels

The [4D](verifications/dim4/) and [5D](verifications/dim5/) branches use radial Chebyshev features and cyclic signs, with nonnegative weights on published near-neighbor count bounds. Their extra objective term is the corresponding weighted count cost.

The dimension-specific external inputs come from Table 2 of [de Laat, de Oliveira Filho and Vallentin, Upper bounds for packings of spheres of several radii](https://arxiv.org/html/1206.2608). Those bounds concern average degree in an epsilon-tangency graph with half-open distance intervals. The new witness replay does not reprove the original spherical SDP inputs. An average count must not be rounded down as if it were an integer count at every center.

### 5. Gaussian-damped cyclic kernel in 48D

The [48D package](verifications/dim48/) uses a separable radial Gaussian and a Gegenbauer-positive angular polynomial. Its global cyclic sign combines a single-sign-change lemma, compact triangle interval checks, and analytic tails. It does not use the finite radial-hat support of the newer 17D–20D family. Changing dimension requires recomputing the angular coefficients, damping threshold, Fourier pair, coupling, and objective; the proof constants do not transfer automatically.

## Exact origin saturation: a common postprocessor

For an isolated-origin degree-zero block, write

    A = [[a,b^T],[b,B]],    alpha = a - b^T B^-1 b.

If B is positive definite and A is PSD, then alpha>=0 and replacing a by b^T B^-1 b gives the exact PSD factorization

    A - alpha e0 e0^T = W B W^T,    W = [(B^-1 b)^T; I].

Only K(0,0) changes. Every nonzero pair, cyclic, far, Fourier, and exterior constraint is unchanged. The density decrease is ball_volume * alpha/H0. This is optimal with B and b fixed, not a new multipoint theorem. See the [earlier origin-slack proof](verifications/origin_slack/) and the [cyclic saturation package](verifications/cyclic_origin_slack/).

The saturated block is singular. A verifier that insists on strictly positive pivots for the full block will reject this valid PSD representation; check its exact factorization instead. The 48D rank-one degree-zero representation is already Schur-saturated, so this particular operation cannot improve it.

The discovery programs impose a strict PSD floor before radial scaling. An interior margin can leave appreciable origin slack despite a small nominal floor. For future candidates, perform exact saturation before ranking objectives, and retain both the source and transformed representation.

## Targeted repair of a failed radial tail

The wider-support 17D polynomial witness originally had four forbidden exterior roots. Its uncorrected objective was not a bound. The [repaired 17D package](verifications/dim17/refinement_repaired/) adds, with X=pi|x|^2,

    g(x) = eta (c-X) exp(-X/4),
    ghat(xi) = eta 2^n (c-2n+16Y) exp(-4Y),    Y=pi|xi|^2.

The exact parameters are c=667/100 and eta=13/10^15. Because c<pi*s, the correction is negative at every permitted nonzero packing separation and preserves the already-verified pair inequality. The matrices and triangle functions are unchanged.

The Fourier correction can be negative at low frequencies, so it cannot simply be assumed positive. Exact polynomial checks control that interval, and its sign is nonnegative beyond the explicit cutoff. Rational root counts and interval bounds cover every physical positive island and the complementary infinite tail. The new Fourier origin is H0'=H(0)-eta*2^n*(2n-c), and the corrected numerator includes eta*c. Subsequent origin saturation uses alpha/H0', not alpha/H(0).

This repair has been independently checked for this witness. It does not repair an invalid PSD block, a positive cyclic sum, or an arbitrary failure of Fourier positivity. A polynomial-only verifier cannot validate the complete mixed-Gaussian function.

## What transfers, and what does not

- Reusable components include exact matrix factorization, Fourier recurrences, Sturm sequences, Bernstein conversion and subdivision, radial Taylor bounds, and the registry schema. Numerical coefficients remain dimension-specific: changing n changes both Fourier and Gegenbauer parameters and the volume multiplier.
- The old hat span and the solid-harmonic hat span differ. Multiplication by u^(ell/2) is not a fixed congruence within the old finite radial span. Their finite models need not dominate one another. Thus the noncompetitive 21D solid-cyclic output does not show that cyclic relaxation itself is weaker.
- Per-harmonic fixed node rescaling within the solid family is a coordinate change in exact arithmetic. However, its numerical strict-PSD floor is coordinate-dependent. Conditioning changes must be distinguished from changes in the underlying function family.
- The generic solid verifier currently asserts 3<=n<=48 and s<T<2s. These are implementation limits, not impossibility theorems.
- The narrow 17D–20D audits pass enlarged rectangular domains without geometry exclusions. The wider repaired 17D audit instead uses three exact exclusions of impossible far-edge triangles. Its proof records the rational inequality establishing each exclusion. Neither audit licenses unverified exclusions elsewhere.
- Exact endpoints matter: the 18D truncation square is 26999999999999997/10^16; the wider 17D value is 29749999999999996/10^16. Replacing them by nearby simple rationals changes the witness. Future parameter grids should start from intended rational strings.
- Kissing configurations and lower packing constructions are not upper-density certificates. They may suggest useful geometry or a separately proved local inequality, but they cannot replace the global packing argument.

## Finite comparisons already performed

| Experiment | Outcome | What follows |
|---|---|---|
| Narrow cyclic 17D–20D | Exact witnesses and independent domain audits pass | Valid source bounds, refined as above |
| Wider 17D without correction | Exact exterior sign fails | Rejected original witness; numerical objective is not a bound |
| Wider 17D with explicit repair and saturation | Independent Fourier, pair, compact, far, saturation and objective checks pass | Refined 17D certificate |
| 17D same-family pointwise controls | Numerical candidates only | No certified optimality gap or isolated attribution of gains |
| 21D solid-cyclic candidate | Numerical objective above the existing registered result | No improvement demonstrated; no infeasibility theorem |
| 22D/23D origin-only screen | Exact degree-zero saturation gives unverified objectives about 0.003914245209 and 0.002901867100 | Still noncompetitive; no global proof claim |
| 25D wider candidate | Retained numerical result; exact degree-zero annular block has a negative elimination pivot | Rejected PSD gate; origin-saturation lemma does not apply |

The last three rows are research diagnostics, not registered density bounds. Files missing from a retained snapshot are not evidence of a running experiment or successful completion.

## Next experiments: bounded tests and acceptance gates

Everything in this section is proposed research, not a verified result.

### Priority 1: certificate-aware radial optimization

Use one promising kernel with verified PSD and triangle signs, freeze it, and optimize only its radial Fourier pair. This isolates radial failures and is cheaper than another full kernel search. Compare an ordinary objective solve with a second solve that caps the objective and maximizes a verification margin.

Near-optimal objective caps and rational Gram-factor reconstruction are established practice; see Section 5.3 of [the several-radii paper](https://arxiv.org/html/1206.2608). The current normalized 10^-10 sign margin is close to solver feasibility tolerances, so a tiny numerical residual can become an actual forbidden root after exact export.

**Gate:** exact PSD/Fourier checks and complete one-variable signs must pass before expensive triangle verification. A changed Gaussian mixture must verify its full transformed function and normalization, not just its polynomial part.

### Priority 2: a six-solve local parameter batch

Take fixed 17D and 18D baselines. For each dimension, run three one-change variants: separation scale s multiplied by 19/20; s multiplied by 21/20; and one additional radial knot placed near the weakest certified margin. Keep all other design choices fixed, preserving the support ratio when changing s. These are six discovery solves plus baseline replays, not a broad dimension sweep. Use exact rational parameter values.

The convenient choice s=n/8 is not a theorem. At fixed normalized separation, changing s tunes Gaussian decay. Knot insertion can represent the previous piecewise-linear function exactly; verify that embedding before interpreting numerical changes. Zero-padding higher-degree blocks similarly nests the mathematical family but may conflict with imposed positive-definite floors.

**Gate:** compare only fully verified, saturated bounds against the corresponding baseline; also record solve time and verification cost. Stop the batch if exact failures show the same conditioning issue rather than a new geometric opportunity.

### Priority 3: isolate basis and cyclic effects

Choose one fixed dimension, support, radial span, and angular degree. Solve the individual-sign and cyclic-sign versions under identical export and certification rules. Both feasible families and all retained comparison witnesses must be clearly identified.

**Gate:** two valid witnesses can compare achieved bounds, but their difference is not the exact difference between the formulation optima. Claiming an optimum gap would require further optimization certificates.

### Priority 4: hybrid features or 48D transfer

A direct-sum kernel combining old hats and solid-harmonic features can contain both families when zero blocks are allowed. Its larger matrices and higher-degree denominator-cleared triangle polynomials may make certification expensive. Prototype at low angular degree before expanding.

For the 48D branch, test a small radial-feature or damping extension with joint reoptimization. A convex combination of complete normalized certificates cannot beat its better endpoint; additional freedom and fresh global checks are needed.

**Gate:** first recover the old certificate exactly inside the enlarged model, then require a competitive complete proof. No cross-dimensional or worldwide-priority conclusion follows from a successful numerical run.

## Research and release discipline

Keep discovery, proof, and registration separate:

1. Retain atomic numerical checkpoints, exact intended parameters, seeds, dependency versions, and residuals.
2. Export immutable exact witnesses. Record analytic inputs and every computational proof gate. Store actual failures such as a positive triangle, forbidden root interval, or negative pivot. An exhausted subdivision budget is inconclusive, not a counterexample.
3. Run the independent audit, recompute the outward objective, update [RESULTS.json](RESULTS.json), and run the complete repository replay. Keep superseded witnesses and their provenance.

The 17D–20D comparisons are pinned to [Afkhami-Jeddi et al.](https://arxiv.org/abs/2006.02560v2) and [Cohn's general-packing table](https://cohn.mit.edu/sphere-packing/). [The newer 2512.04829 paper](https://arxiv.org/html/2512.04829v2) reports the ranges 4–7 and 9–16, not 17–20. Indexed freshness and unreplayed external comparison data do not establish worldwide priority. Percentages describe reductions of upper bounds, not improvements in the unknown optimum.

## Exact transfer to dimension 47

A further verified route restricts the Gaussian-damped 48D physical function to dimension 47. Fourier positivity follows from its exact one-coordinate marginal, independently reconstructed coefficient by coefficient. The same angular polynomial has nonnegative dimension-47 Gegenbauer coefficients and a larger constant term. An exact rank-one origin Gram saturation lowers its origin cost while leaving every required nonorigin and cross value unchanged. All inherited continuum and tail proofs remain required.

This gives density <6.16879891e-7, about 4.96% below the pinned published comparator. The same direct-transfer screen is noncompetitive in dimensions 25–46; that screen is not an impossibility theorem. Numerical Gaussian and integration perturbation screens are research only. See [the full transfer proof and reproducible code](verifications/dim47/PROOF47.md).

## Subsequent verified branches: 6D, 45D and 46D

The 6D count-enhanced Chebyshev branch gives <0.410189987 after exact Gaussian guards and origin Schur saturation. It retains external published near-neighbor count hypotheses, which are not newly reconstructed. All new continuous-domain conditions are checked. This small improvement is distinct from the noncompetitive low-dimensional solid-harmonic trials.

For 46D, a new degree-111 radial Fourier pair combined with the inherited angular polynomial gives <9.14882847e-7. For 45D, changing both the angular polynomial and radial pair gives <1.311423575e-6. Both require fresh global radial signs and coupled pair bounds; 45D additionally requires a fresh angular/cyclic proof. Exact dimension-specific origin Gram saturation remains essential. These results supersede the earlier noncompetitive direct-transfer screen only for these newly constructed witnesses, not by reclassifying the old experiments.

See the complete [6D](verifications/dim6/), [45D](verifications/dim45/) and [46D](verifications/dim46/) packages, and the pinned source comparisons in [COMPARISONS.md](COMPARISONS.md).

## Six-dimensional continuation refinement

A separately verified continuation witness improves 6D to <0.409135399. Unlike the earlier repaired witness, both Gaussian guard weights vanish; exact Sturm checks establish global polynomial signs directly. The search mixed a warm-start continuation endpoint with a near-feasible joint radial/kernel result, but no numerical convex-feasibility assumption is used in the proof. The final exact matrices, Fourier pair, origin Schur operation, count costs and every continuous domain are checked anew. This preserves the external count-theorem assumptions rather than claiming a count-free improvement.

## Dimension 44 and the small six-dimensional follow-on

The new 44D Gaussian-damped cyclic witness changes both angular and radial data, with exact minimum squared separation 23/4 and radial decay parameter 51/23. It gives <1.842868599e-6 after fresh angular, cyclic-domain, global radial-sign, coupled pair and normalization proofs. No old 48D coefficients or signs are substituted for these checks.

The separately audited small 6D follow-on gives <0.409128533. Its final coefficients pass global polynomial signs without Gaussian repair; the existing count assumptions, Schur operation and full-domain conditions remain required. Both previous 6D certificates are preserved.

## Guarded 5D, count-enhanced 7D and tetrahedral research

The guarded 5D result <0.510129239 needs an explicit exact PSD kernel correction: the unguarded numerical checkpoint is invalid. This correction is separate from Gaussian guards. The 7D result <0.320584772 uses the count-enhanced cyclic family with full exact/interval verification. Both retain their cited external count hypotheses.

The research-only [22D tetrahedral package](research/tetrahedral_aggregate/) proves 3T <= BE for squared shell widths q=21/20 and q=11/10, with strict B bounds 3341.14916082 and 5086.98083741. The latter is about 36.75% below the stored affine certificate, a comparison of valid witnesses rather than proved optima. Exact endpoint contradictions separate both new polynomials from the affine-envelope family; a harmonic reduction theorem explains why adding unweighted higher normal harmonics does not fix that restriction. These auxiliary four-point results do not establish a new density bound. Density activation material remains isolated as negative/provisional research and is never a proof input or registry entry.

## Five-dimensional long-edge lens extension

The lens witness <0.509996147 expands the admissible kernel family rather than merely retuning coefficients in the old cone. It permits a small positive far contribution on triangles with two short sides and one uniquely eligible long base. Elementary midpoint-ball centroid bounds cap common neighbors of that base. The physical tail supplies a nonnegative step allowance, so exact incidence counting pays every permitted triangle contribution through the pair term. Both pair and triangle contributions carry the same factor two.

The proof works on simultaneous-translation orbits for sufficiently large periods and then uses periodic approximation, so motif images are not omitted. Step endpoints receive the minimum adjacent allowance; outside support the allowance is zero. Full independent checks establish global Fourier signs, PSD, pair inequalities, compact cyclic signs, every long-base budget bin and the newly relaxed far domain. An exact positive far witness shows why the old theorem would not suffice. The elementary lens caps add no external spherical-code assumptions; published near-neighbor counts remain as before.

The stronger spherical-lens continuation gives <0.509953498. It improves the same long-edge accounting through conservative transverse projection, internally replayed exact spherical LP caps and Rankin/centroid bounds. This adds exact cap verification and changes the tail-budget support; it does not treat numerical spherical-code estimates as assumptions. The earlier centroid certificate and the positive old-far witness remain available to make the theorem change explicit.

## Gaussian cyclic continuation in 41D–43D

The new [41D](verifications/dim41/), [42D](verifications/dim42/) and [43D](verifications/dim43/) certificates continue the Gaussian cyclic family used in 44D. Restriction alone, even with origin saturation, did not improve the pinned upper bounds in the bounded 22–43D transfer screen. The successful change was to optimize a new dimension-specific angular polynomial together with a degree-31 Laguerre correction to the earlier degree-111 radial function.

The angular polynomial has nonnegative normalized Gegenbauer coefficients, a certified negative cap and a single subsequent sign change. Complete cyclic-domain subdivision, plus analytic scale and tail arguments, proves the required triangle sums. A rank-one origin Gram block supplies the optimal origin cost for that fixed angular/radial form. Exact Fourier reconstruction, global ray signs and complete radial coupling finish the proof. No external neighbor-count input is used.

The first sampled candidates were not valid certificates. Exact angular separation and root-guided radial cuts were essential. Discovery scripts and frozen checkpoints are separate from the self-contained proof packages at [41D discovery](research/dim41_gaussian_discovery/), [42D discovery](research/dim42_gaussian_discovery/) and [43D discovery](research/dim43_gaussian_discovery/). Their numerical feasibility is never used as proof.

Gap comparisons now use the newer lower constructions in arXiv:2607.20359v4, rather than the older indexed table. See [comparison conventions](COMPARISONS.md) for the formulas and pinned upper references. A finite successful continuation does not establish a universally improving method in all dimensions.
