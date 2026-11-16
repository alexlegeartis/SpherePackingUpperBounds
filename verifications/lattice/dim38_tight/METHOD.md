# 38-dimensional lattice cyclic-collinear hybrid

Safe upper density: 0.000013058679142.
Frozen certificate: screen_n38_collinear_fixed_tight_certificate.json.
SHA256: 34888f6e673b2369bf69cd506e6fab0c3cf8612a03c3617ff85ca0d1678a9eda.

Suggested table label: **Cyclic PSD kernel with forced collinear triples + mod-2 class profiles**.
The proof additionally uses cross-class pair-sum angular inequalities and nested conditioned Gram-moment class capacities.

## Complete proof obligations
1. Published dimension38 angular polynomial Q and profile h(t)=t^96 exp(-120(t-1)): origin-augmented positive semidefinite kernel and cyclic triangle inequality, including noncompact tail and all dilations. Both published source gates freshly replayed in lattice_cyclic38_replay; these are essential dependencies, not a previous density-bound assumption.
2. Lattice cyclic summation and collinear orbit identity. The ratios {-1,2,1/2} form a disjoint cyclic-invariant ordered-pair subset. Its contribution is sum_v[K(v,-v)+2K(v,2v)]. Thus sum_v(D(v)+R(v)) >= -C. Absolute convergence follows from Gaussian decay. See independent CYCLIC_LATTICE_BRIDGE.md and COLLINEAR_RESIDUAL.md.
3. Rational normalized Laguerre reconstruction, q nonnegative on [0,infinity), p nonpositive after the envelope cutoff. Exact global Bernstein proofs saved separately and independently replayed.
4. Mixed spatial envelope on EVERY bin: f(t)+C/B[(1+Q(-1)/B)h(t)^2-2h(t)+2h(t)h(4t)] <= g_bin. Positive cross term is retained. Frozen producer uses outward Arb Taylor16 plus full Leibniz remainder; independent verifier uses a separate centered Taylor/differential-operator implementation. The combined tail is nonpositive.
5. Exact pair-sum angular inequalities: positive normalized Gegenbauer coefficients, continuous sign/majorant checks, valid projected projective cap polynomials, whole-bin source indicators and all target bins. Four active cuts are embedded.
6. Per-class cumulative bounds A/B, unconditional nested weighted rank constraints, narrow-short conditioned rank constraints, and all-short constraints with slope 21/20 times base slope. All-short theorem independently audited for dimension38 and refined cutoffs; retain exact endpoint hypotheses. There is NO density-index exclusion premise.
7. Exact nonnegative rational LP duals dominate every long-bin envelope coefficient for all 193 short-vector or no-short-vector class branches. Their maxima, plus angular tangent costs, bound the lattice sum.
8. Poisson normalization and exact even-dimensional density comparison. Objective is (p(0)+C+2m*kappa+2*sum nu*N0)/q(0). C is not doubled; opposite-pair factor2 applies only to class sums. m=2^38-1.

The prior live comparator 0.0000145182018536 is used solely for the improvement percentage, not in the proof. Do not add an additional forced-line/perfectness credit to this certificate without proving it does not double count the already included collinear orbit.
