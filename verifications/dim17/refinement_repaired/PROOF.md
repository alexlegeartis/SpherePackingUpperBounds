# A certified slow-Gaussian repair of physical-tail sign defects

## Result and status

The wider-support 17D cyclic witness originally failed the exterior physical
sign test. A separate, analytically controlled Gaussian correction repairs it
without changing any matrix. The complete primary replay now proves

Delta_17 < 0.018250915.

Its density enclosure is centered at
0.01825091405434677552501772408911594661720864933457687712605626978655699664777,
with radius below 5.68e-78. Independent mixed-Gaussian and Schur review passed; see REVIEW.md. This is
not a claim of a new general sphere-packing theorem or worldwide priority.
The underlying cyclic density argument is unchanged from [the original cyclic proof](../PROOF.md).

The source is certificate.json, and the correct verifier is
verify_repaired.py with slow_gaussian.py. The original polynomial-only verifier
must not be used to interpret the corrected Fourier pair or objective.

## General repair lemma

Write X=pi|x|^2, Y=pi|xi|^2 and a=n/2. Suppose P(X)e^(-X) has Fourier
transform H(Y)e^(-Y). Assume H>0 globally, the PSD kernel and its cyclic
inequalities are already valid, and the base coupled pair inequality holds.
Only the physical-tail sign is permitted to fail.

Choose rational c with 0<c<pi*r^2 and eta>0. Add the correction

g(x)=eta*(c-X)*exp(-X/4).

Differentiating the Gaussian Fourier identity in its width gives exactly

ghat(xi)=eta*2^n*(c-2n+16Y)*exp(-4Y).

This identity holds in odd dimensions as well: (1/4)^(-n/2)=2^n is rational.
Since c<pi*r^2, g is strictly negative at every permitted nonzero packing
separation. Thus adding g preserves the coupled pair inequality; the matrices,
PSD property, compact triangle inequalities and far-edge inequalities are
unchanged before normalization.

Put D=2^n*(2n-c)>0 and y0=(2n-c)/16. For Y>=y0, ghat>=0. On [0,y0], after
dividing by exp(-Y), the correction is bounded below by -eta*D. Consequently
it suffices to prove the exact polynomial inequality H(Y)>eta*D on [0,y0].
This can be checked by rational Sturm root counts and endpoint signs. The new
Fourier origin is N'=H(0)-eta*D, which must be positive.

For X>=pi*R^2, the physical sign is equivalent to

P(X) <= eta*(X-c)*exp(3X/4).

Where P is already negative this follows immediately. Cover every remaining
positive island by finitely many rational intervals [l,h], with l>c. Exact
Sturm checks establish negativity on all complementary intervals and the
infinite final ray. On each repair interval, a certified upper bound U for P
suffices if

U < eta*(l-c)*exp(3l/4).

The right-hand side is a lower bound throughout [l,h]. This proves the entire
infinite physical tail without replacing global signs by sampling.

Finally divide both f+g and the unchanged unnormalized kernel by N'. The
corrected density objective before the ball-volume multiplier is

[P(0)+eta*c+A_0[0,0]] / [H(0)-eta*D].

An optional subsequent exact origin Schur subtraction alpha changes only the
numerator to P(0)+eta*c+A_0[0,0]-alpha. The decrement must use N', not H(0).

## Exact 17D instance

Use the seven-node, angular-degree-16, radial-degree-31 wide witness. Minimum
separation squared is 17/8, and the upper squared cutoff is the exact stored
rational 2.9749999999999996. Do not silently replace it by 119/40.

Take beta=1/4, c=667/100, eta=13/10^15. Then

eta*D = 0.00000004656857088,
y0 = 2733/1600.

The exact polynomial H-eta*D is positive throughout [0,y0]. Physical repair
intervals are [2113/50,423/10] and [6109/100,1223/20]. Arb bounds respectively
show

- P <= 15.06585467..., below correction allowance 26.92991261...;
- P <= 52662470.398..., below correction allowance 55974382.904... .

Exact Sturm tests prove P<0 on every complementary interval and the infinite
last ray. The old tail had four forbidden real roots, so it was invalid until
this correction was proved.

The full replay also verifies exact rational PSD, the original Fourier identity,
the coupled pair inequality in 42 boxes, all 56 compact radial-cell triples and
all 21 far radial-cell pairs. Those geometric checks are not inferred from the
narrower earlier witness. See full17_repaired.log.

## Reproduce

python verify_repaired.py certificate.json

Never use Python -O. Verification dependencies remain NumPy and python-flint.
No optimizer or network is needed. The proof uses exact rational arithmetic,
Arb interval bounds and the explicitly documented outward-rounded Bernstein
subdivision. It is not a formal proof-assistant verification.

## Why this helps and its limits

Numerical optimization can leave tiny normalized violations that become
positive islands in a large polynomial at large X. Their Gaussian-weighted
values can still be extremely small. A slower-decaying negative physical
correction can dominate those islands with very little density cost, while its
Fourier negative part is confined to a short low-frequency interval.

This does not repair an invalid PSD kernel, a positive cyclic triangle sum,
a failed coupled pair inequality, or arbitrary Fourier negativity. Those remain
separate mandatory proof gates. It also does not justify disregarding any
positive tail island, however small. The correction and changed normalization
must be explicitly part of the exact witness and verified.

## Registered origin-saturated conclusion

The source above remains positive definite and proves the unsaturated bound <0.018250915. The registered <0.018204096 is obtained by the exact Schur derivation in REVIEW.md and repair_checks.py: subtract alpha e0 e0^T from the origin block, reconstruct the resulting PSD matrix as W B W^T, and decrease the numerator by alpha while keeping the corrected mixed-Gaussian Fourier normalization. The independent rational Machin bound explicitly proves the saturated strict decimal. All off-origin constraints and the three exact geometric exclusions remain unchanged. The unchanged source certificate plus this derived transformation specifies the registered witness; no singular matrix is passed to the strict-PD source check.
