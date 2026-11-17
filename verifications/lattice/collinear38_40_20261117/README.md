# Sharper primitive-line decay and optimized PSD origin, dimensions 38–40

These are small certified refinements, not a new large-margin breakthrough or an exhaustive worldwide-priority claim. All external comparisons retain Zubrilina (2026), [Table 1](https://arxiv.org/html/2610.35830v1#S4.T1). Earlier repository bounds are separately identified as internal history.

| Dimension | Certified density upper bound |
|---:|---:|
| 38 | <0.000012713806752 |
| 39 | <0.000008757804412 |
| 40 | <0.000006157949542 |

## Proof
Use precisely the frozen Q42 polynomial and predecessor radial/profile certificates in [rank_two_20261117](../rank_two_20261117/). Their complete rank-two sign theorem, nonnegative target-dimensional Gegenbauer expansions, Fourier positivity, physical majorants and mod-2 budgets are unchanged.

Let h(t)=t^96 exp[-120(t−1)] and h_k=h(k²t), t≥1. Then
h_(k+1)/h_k=((k+1)/k)^192 exp[-120(2k+1)t]≤exp(192−360)=exp(−168)<ρ=10^−60.
The first inequality uses log 2<1, k≥1 and t≥1. The final inequality is certified by an exact positive Taylor lower bound for exp(168). Also 0<h_k≤1 because h decreases on [1,∞).

For i≥1 and j>i, h_i h_j≤h_1 h_2 ρ^(i+j−3). Summing the geometric series gives Σ_(i<j) h_i h_j≤h_1 h_2/Γ, Γ=(1−ρ)(1−ρ²). This is the same predecessor full-line domination argument, but with a much sharper rigorously justified ratio.

Set c=C/B², u0=c q_old, v0=c R from the saved predecessor. Use u=cBΓ/2, v=v0+u−u0. Exact arithmetic verifies u≥u0>0. The difference between new and old complete line corrections is
4(u−u0)Σ(h_k²−h_k)+8uΣ_(i<j)h_i h_j−4cBΣ h_k h_(2k)≤0,
because h_k²≤h_k, 2u=cBΓ and Σ h_k h_(2k)≥h_1h_2. Absolute Gaussian convergence justifies the infinite sums.

Normalize the new angular polynomial to constant Gegenbauer coefficient 1 and put q=Q(1). The origin block is PSD with origin scalar q v²/u. Thus the existing Poisson/PSD/profile proof applies without any new sampling assumption or new extremality premise. Exact rational objective-ratio conversion from the predecessor's certified safe upper bound is conservative, including odd dimension 39. Every displayed new decimal is strictly above that exact rational upper bound.

## Replay and audit scope
Install the versions listed in the unchanged predecessor requirements. Run:

```sh
python verifications/lattice/collinear38_40_20261117/verify.py
python verifications/lattice/collinear38_40_20261117/verify.py --full
```

The quick verifier independently reconstructs target Gegenbauer expansions, checks frozen hashes and all transfer inequalities, and regenerates RESULTS.json with exact rational parameters. The full command additionally replays all three complete predecessor proofs, including the Q42 interval gate and universal tail bounds. The inherited numerical interval implementation is rerun, not independently reimplemented for this refinement. Prior proof packages and their historical registries are retained unchanged; the live LATTICE_RESULTS.json points here.

Fresh complete selected predecessor replay passed for all three dimensions; see `predecessors_full.log`. Q42 closed all 16,093 interval nodes with no pending boxes. All radial and profile checks passed. Exact transfer verification also passed independently.
