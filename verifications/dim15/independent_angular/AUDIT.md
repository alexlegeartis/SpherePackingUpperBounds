# Independent 15D curved-profile angular and theorem audit

**PASS.** This audit covers the angular positive-semidefinite kernel, its exact sign structure, the complete cyclic geometric domain, the all-scale extension, and the finite-sum theorem linkage. It does not audit any radial polynomial or assert a density bound on its own.

Pinned mathematical input:

- `inputs/angular_exact.json`
- SHA256 `ebeb0fd1b0532445cf274f01f6f01f1996f608fc9794eb09e3cb43f3932f6354`
- Dimension15, normalized Gegenbauer constant coefficient1, all Gegenbauer coefficients nonnegative.
- Profile h(u)=exp[-12(u/s-1)-16(u/s-1)^2], for squared radius u>=s, with s=15/8.

The same exact Q appears in `../dim15_curved_ready_20261111`. No file in that package or the discovery folder was edited. The changing radial discovery inputs were not read for this audit.

## Independent computational checks

`independent_check.py` reads only the pinned exact angular input. It does not import the primary verifier or regenerate its coefficients from discovery floats.

1. **Exact reconstruction.** It reconstructs the normalized Gegenbauer expansion using the closed coefficient formula, independently of the primary three-term recurrence. The rational monomial polynomial, coefficient signs, constant harmonic coefficient and B=Q(1) agree exactly. Trailing zero harmonics are allowed; this Q has actual polynomial degree12.
2. **Exact cap.** Exact rational Bernstein coefficients on20 subintervals prove Q(t)+1/20000<0 for every t in[-1,1/2], without using the primary Sturm implementation. The exact leaves are saved in `cap_partition.json`.
3. **One sign transition.** Complete320-bit Arb isolation accounts for every complex root with multiplicity. Exactly one possible real root occurs in(1/2,1), with endpoint signs Q(1/2)<0<Q(1). Therefore there is exactly one real transition, at approximately0.50000049267291159545. The exact cap already rules out roots on[-1,1/2].
4. **Manual derivative audit.** SymPy verifies all eight angle/exponent partial derivatives used by the primary code, as exact symbolic identities.
5. **Distinct cyclic interval proof.** A separate automatic-differentiation implementation and a different rational box subdivision rule prove the entire compact cyclic domain at320 bits. It accepts10,744 rectangles and discards1,041 rectangles solely by exact sorting exclusion. All23,569 tested nodes finish, maximum binary subdivision depth29. Accepted and excluded rectangles have exact total area9/2, equal to the initial rectangle. The full partition and rigorous upper enclosures are saved in `interval_partition.json`.
6. **Analytic tails.** Both unbounded-domain inequalities pass with wide margins. The two positive upper bounds are approximately1.26414e-1638 and1.25990e-12, versus the exact cap margin1/20000.

The original verifier was also copied and replayed in `primary_replay/`. It reproduces the frozen mathematical angular input byte-for-byte and passes with5,752 accepted rectangles,624 sorting exclusions and per-axis depth15. That supplemental replay is not needed by the independent entry point.

## Geometry and derivatives inspected

Write the normalized sorted side lengths as1, U, W, with U>=1 and W>=U. If T is the cosine opposite W, then

W^2=1+U^2-2TU, -1<=T<=1/2, and 2TU<=1.

These parameters cover every triangle whose shortest squared side equals s, including degenerate collinear limits. Its other cosines are c1=(1-UT)/W and c2=(U-T)/W. Sorting gives T<=c1<=c2.

The denominator is uniformly positive on the entire rectangular domain, including its unsorted part:

W^2=1+U(U-1)+U(1-2T)>=1.

This identity justifies the independent checker's explicit refinement of the square-root input to its analytic lower bound1. Arb ball multiplication can otherwise lose one-sided positivity through dependency; it does not invalidate the bound. Reciprocals are taken only after a strictly positive denominator is established. Derivatives are enclosed by interval automatic differentiation, not copied from the primary manual formulas.

Set d=1-2TU and e=U^2-2TU. The two logarithmic profile ratios are

E1=-12d-16d(2U^2-2+d), E2=-12e-16e^2.

The symbolic audit verifies the four cosine partials and these four exponent partials:

- E1_T=-(12+32e)d_T.
- E1_U=-(12+32e)d_U-64dU.
- E2_T=-(12+32e)e_T.
- E2_U=-(12+32e)e_U.

The evaluated expression is Q(T)+exp(E1)Q(c1)+exp(E2)Q(c2). The centered polynomial enclosures clip cosine ranges to[-1,1]. This is safe for every point of the full rectangular box, not only its sorted subset: T specifies a real angle, U>=1, and W is the corresponding third side. Thus the mean-value enclosure remains valid even when a box center violates sorting. The proof never treats a finite sample or a gradient approximation as an acceptance certificate.

Every discarded box has2*T_lower*U_lower>1. This can happen only with positive T_lower, so all points in it violate sorting. Other boxes are accepted only by a strictly negative rigorous direct enclosure or a rigorous gradient mean-value upper bound. Rational binary subdivisions cover their parents without gaps; no unresolved box remains.

## Unbounded U>=4 tails

Sorting implies W>=U and e=W^2-1>=15. Thus exp(E2)<=exp(-3780).

If T>=-1/4, the numerator1-UT is nonnegative because2TU<=1. Hence c1<=(1-UT)/U=1/U-T<=1/2, and its Q term is nonpositive. It is sufficient that B*exp(-3780)<1/20000.

If T<=-1/4, then d>=1+U/2>=3 and2U^2-2+d>=0, so exp(E1)<=exp(-36). It is sufficient that B*(exp(-36)+exp(-3780))<1/20000.

Both checks are rigorous Arb inequalities. The bound |Q(t)|<=B follows from its nonnegative normalized Gegenbauer expansion, equivalently from the two-point PSD kernel inequality. No geometric tail is omitted.

## All-scale extension

For sorted squared sides x<=y<=z, all at least s, let q0,q1,q2 be Q at the opposite angles in the order z,y,x. The largest angle has cosine at most1/2, so q0<0. The one-sign-transition property and angle ordering allow only the secondary sign cases(--),(-+),(++), with zeros included.

After division by the positive h(x)h(y), the cyclic sum is

q0+r1*q1+r2*q2, where r1=h(z)/h(y), r2=h(z)/h(x).

For v>=u>=1 and common squared-side dilation lambda>=1, the logarithmic derivative of h(s*lambda*v)/h(s*lambda*u) is

-(v-u)*(12-32+32*lambda*(v+u))<=0.

The bracket is at least44. Consequently r1 and r2/r1 both decrease under dilation. If both secondary q values are nonnegative, their weighted terms decrease. If both are nonpositive, the total is immediately negative. In the mixed case q1<=0<=q2, write the secondary sum as r1*(q1+(r2/r1)*q2). Its bracket decreases; while positive the product of its two nonnegative decreasing factors decreases, and after it crosses zero the secondary sum remains nonpositive. Therefore the cyclic inequality at shortest squared side s extends to every triangle with sides at least sqrt(s).

This argument handles equal sides and zero q values. It does not assume that negative summands individually decrease, and it makes no false termwise-monotonicity inference.

## PSD origin extension and finite-sum linkage

Conditionally on M>=0, define the off-origin kernel as

K(x,y)=(M/B)h(|x|^2)h(|y|^2)Q(cos angle),
K(0,x)=-M h(|x|^2), K(0,0)=MB.

The degree-zero feature is the rank-one Gram form with origin feature B and nonzero feature -h; all higher harmonics are PSD because their coefficients are nonnegative. Nonzero radii below sqrt(s) can be assigned zero features without changing any separated-star condition. The pair diagonal is K(x,x)=M h(|x|^2)^2, giving exactly the required pair inequality f(x)<=M(2h-h^2).

For a finite set, sum each rooted-star Gram form. Ordered pair contributions are2K(0,x)+K(x,x); each unordered triangle contributes twice its cyclic sum. The latter is nonpositive. Applying the pair inequality therefore gives sum_{i,j} f(x_i-x_j)<=N(f(0)+MB), with no missing factor of two. For periodic packings, Gaussian-polynomial f is Schwartz and its nonnegative Fourier transform gives the usual Poisson lower bound. The curved h decays as exp(-constant*|x|^4), so all separated-star sums converge absolutely and boundary terms vanish. This yields the stated objective and volume factor once a separate radial audit supplies Fourier positivity, H(0)>0, pair inequalities and the exact objective. This audit does not supply those radial premises.

The optional inverse-power obstruction in the source theorem is also correct: on collinear side lengths1,r,r+1, a decaying profile with h((r+1)^2)/h(r^2) tending to1 has limiting normalized cyclic sum Q(-1)+Q(1)=2*sum_even c_l>=2c0>0. It is not needed for the accepted bound.

## Reproduction and scope

Run `verify_all.py` with Python3, python-flint and SymPy. In this workspace:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 research/density_solid_cyclic_20261111/.venv/bin/python research/density_curved_profile15_angular_independent_20261111/verify_all.py
```

The independent interval pass takes approximately6 seconds here. Reproduction regenerates result/partition files; copy the folder first if preserving their original hashes. The input hash is enforced by the checker. `SHA256SUMS` records the frozen audit files.

No radial candidate, finite-grid optimization claim, or density rounding is accepted by this audit. A density conclusion requires the separate radial certificate with exactly this Q and B.

## Assertion-mode packaging correction

Both standalone entry points and the supplemental primary replay explicitly reject assertion-disabled execution before loading proof inputs. `test_optimization_guards.py` verifies rejection under both `-O` and `-OO` for all three scripts. The normal independent full replay and supplemental primary replay pass after this guard-only correction. The pinned angular mathematical input is unchanged.
