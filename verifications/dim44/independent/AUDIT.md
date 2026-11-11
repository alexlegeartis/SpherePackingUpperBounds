# Independent review of the frozen 44D cyclic density certificate

## Outcome

**PASS: the unrestricted packing density satisfies Delta_44 < 0.000001842868599.**

The independently reconstructed interval is centered at

1.8428685980689728749847153538789270688874740498325860280633221295411415419001698155391739524361615350977790956266403697981038487072882265827373653111156e-6,

with radius below 1.5e-158 at 512-bit precision. Exact rational rounding, independently of that displayed interval, proves the stated strict decimal. Relative to the pinned comparator 0.0000019017703144, the conservative decimal is approximately 3.0972% smaller. This review does not establish worldwide priority.

## Frozen inputs

The reviewed package is `research/dim44_gaussian_ready_20261111/`. Source SHA256 values:

- `certificate.json`: `4ca48c456c914770630100c4f3865967a5a6e039c7a29d78ebadaca2b77c5c55`.
- `continuum_three_point/cyclic_angular40_exact.json`: `2515f35922ff1c78dcc913b7f98ebadf10c924b1c795deb9fce0c80aab0c0b7e`.
- `verify.py`: `f20b07d04ea56cf2762faf1ae62cd259ccc84586262b493438c56680e442bdce`.

These hashes were checked after the matched final export. The complete replay reconstructed the stored objective against this angular polynomial; an earlier source-version issue and earlier rejected numerical candidates are not proof inputs.

The exact parameters are dimension 44, minimum squared separation s=23/4, physical tail threshold 7, and kappa=51/23. The radial polynomials P,H have degree 111. The angular polynomial has degree 19, constant normalized 44D Gegenbauer coefficient exactly 1, and B=Q(1)>0. M is the positive rational in the witness.

## Verification performed

### Complete isolated replay and source review

The complete verifier was rerun in an isolated copy, with assertions enabled, and passed. The review checked the code and mathematical coverage of:

- Exact angular reconstruction, all nonnegative 44D harmonic coefficients, and their sum B.
- Strict angular signs on [-1,1/2] and [-1,5001/10000], plus the unique remaining root before 1.
- The actual cyclic expression on the full ordered compact triangle domain: 6,004 certified rectangles, 418 exact exclusions for the incompatible ordering condition, maximum depth 15. Excluded rectangles contain no point of the required domain; boundary equality remains covered.
- The single-sign-change scaling lemma and analytic tails outside the compact domain. The exact inequality 40^30/30!>2B*100000 controls positive tail contributions. The actual Gaussian scale is s*kappa*pi=51*pi/4>40.
- Complete certified degree-111 root isolation for H on the nonnegative ray and P beyond a rational point below 7*pi, with exact endpoint signs. Each root list is checked for complete degree count. These imply both full infinite-ray signs, rather than sampled positivity.
- The source radial coupling proof on the whole interval, with 332 accepted intervals of maximum depth 13.

The global cyclic and root-isolation implementations were reviewed and rerun, not replaced by a wholly clean-room implementation. This scope distinction is intentional.

### Independently implemented checks

`independent_radial.py` implements separate checks on the same exact witness:

1. **Fourier identity.** Decompose P exactly in the normalized Gaussian-Laguerre eigenbasis for dimension 44, flip coefficients by (-1)^k, and reconstruct H. This differs from the producer's direct monomial transform.
2. **Angular PSD.** Recover the 44D Gegenbauer coefficients by spherical-moment orthogonality, using E[t^(2j)]=(1/2)_j/(22)_j and the exact harmonic dimensions. This differs from triangular coefficient extraction and the explicit coefficient formula. The recovered vector equals the stored nonnegative vector exactly, and its constant coefficient is 1.
3. **Origin normalization.** Check the rank-one degree-zero relation (MB)(M/B)=M^2, with M>0. Therefore the kernel values are K(0,0)=MB, K(x,0)=-M h(x), K(x,x)=M h(x)^2.
4. **Radial coupling.** Independently prove

       P(X) exp((kappa-1)X-kappa*pi*s)
         -2M+M exp(-kappa*X+kappa*pi*s) < 0

   on a rational enclosure of [pi*s,7*pi]. A centered second-derivative mean-value bound, with separately constructed derivative polynomials and Arb interval bounds, passes in 17 leaves, maximum depth 6. All subdivision intervals cover the original interval, and no interval is excluded by sampling. Beyond 7*pi, P<0 and 0<h<=1 give the coupling directly.
5. **Final rounding.** A separate rational Machin enclosure for pi proves

       [(P(0)+MB)/H(0)] * (pi_upper*s/4)^22 / 22!
         < 0.000001842868599.

   The exact objective agrees with the stored objective and H(0)>0.

## Mathematical implication and scope

The PSD feature construction uses the dimension-44 addition theorem. A smooth cutoff near the origin defines the kernel globally without changing any required packing separation. The radial function and kernel features are rapidly decreasing, so periodic pair and triangle sums are absolutely convergent.

The complete cyclic signs and pair inequality give the existing cyclic continuum density conclusion of Cohn, de Laat and Salmon, not the literal individual-sign premise of their Theorem 1.1. Periodic approximation extends the result to arbitrary packings. There is no lattice-only restriction or external neighbor-count assumption.

Only the general Gaussian-cyclic scale and tail arguments are inherited from the historical 48D notes. The angular polynomial, dimension-specific PSD decomposition, Fourier pair, minimum radius, coupling interval, and final objective have all been checked for this 44D witness.

This is a computer-assisted audit using exact rational arithmetic and Arb, not proof-assistant formalization or external peer review. It makes no packing-construction or optimality claim.

## Reproduction

Use Python with python-flint installed, without `-O`:

    python verify.py

The entrypoint verifies frozen input hashes, runs the complete reviewed source verifier, and runs the distinct independent checks. The independent script may also receive the frozen package directory as its first argument. Its output is `independent_radial.json`. The isolated source replay is recorded in `complete_replay.log`; the separate implementation has its own `independent_radial.log`.
