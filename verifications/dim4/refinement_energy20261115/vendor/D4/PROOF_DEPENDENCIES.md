# Proof dependencies of the D4 upper bound

This file audits the available proof of upper density `0.629521681653` for
every packing of congruent balls in `R^4`. It is not a proof that the D4
lattice is optimal. The published proof is self-contained in `paper.tex`,
`data/`, and `verify/`; the separate Lean formalization is unreleased and is
not an input.

## Mathematical chain

1. The scalar function and its Fourier-side certificates are stored in
   [`data/v7/scalar.json`](data/v7/scalar.json). Stage 1 checks the coefficient
   identities, tails, pair signs, all 112 cyclic branches, local-energy
   parameters, and the raw objective.
2. The degree-zero Gram block and the origin allocations give a nonnegative
   anchored energy on every admissible finite neighbourhood. The component
   inventory and written arguments are in
   [`verify/stage1_v7/COMPONENTS.md`](verify/stage1_v7/COMPONENTS.md).
3. Stage 2 projects the degree-zero block onto the five cumulative counts,
   checks the twelve affine floors and the three spherical count certificates,
   and enumerates all 165,503 vectors in the stated count domain. Stage 3
   applies the exact finite-count refinement to all 124 group floors. See
   [`stage2_origin/PROOF.md`](verify/stage2_origin/PROOF.md) and
   [`stage3_census/finite_count_gain.md`](verify/stage3_census/finite_count_gain.md).
4. Stage 4 separates the broad graph term. It checks the annular count 151
   and replays the projected S2 certificate in
   [`data/projected_eight.json`](data/projected_eight.json) with
   [`verify_projected_eight.py`](verify/stage4_graph_slack/src/verify_projected_eight.py).
   The exact Legendre/Bernstein replay proves the projected code bound 8, so
   the graph has at most 1208 ordered edges and contributes `3/10^6` of
   anchored energy.
5. Stage 5 replays every leaf of the two kernel-box subdivision trees. Its
   written far-pair lemma and exact symbolic link are in
   [`far_pair_lemma.md`](verify/stage5_kernelbox/proofs/far_pair_lemma.md),
   [`kernel_box_method.md`](verify/stage5_kernelbox/proofs/kernel_box_method.md),
   and [`kernel_link/README.md`](verify/kernel_link/README.md).
6. The written density-transfer argument in `paper.tex` and
   [`conditional_theorem.md`](verify/stage5_kernelbox/proofs/conditional_theorem.md)
   passes the finite anchored inequality to arbitrary packings. The final
   exact comparison is

       U_raw - a_new/8 - 3/8000000 < 629521681653/10^12.

## Exact projected-code certificate

The data file contains a nonnegative combination of normalized Legendre
polynomials with constant coefficient `99999/100000` and exact objective

    89804130521/9999900000 < 9.

`verify_projected_eight.py` reconstructs the Legendre recurrence with rational
arithmetic, maps the closed interval `[-1,223/735]` to `[0,1]`, and certifies
the sign of the negative polynomial by exact Bernstein subdivision, including
the endpoint. The addition formula for normalized Legendre kernels is the
written positive-definiteness input. The classical Schütte--van der Waerden
nine-point theorem remains a citation and alternative route, not a dependency
of the replay.

## Written and cited inputs

The stage programs establish the finite rational and interval certificate
claims. The following parts remain written mathematical arguments or cited
standard facts: the local cap, tent and allocation inequalities; the
positive-definiteness/addition formula for spherical kernels; the interval
lower-envelope soundness argument; and the finite-to-density transfer. Their
scope and conventions are stated in the paper and the stage proof notes.

The independent kernel assembler and exact stage-1/stage-5 link are additional
checks, not formal semantics proofs of Python. Hash manifests establish input
identity, not mathematical correctness.

## Formalization boundary

The Lean formalization is unreleased work and is not shipped with this package.
See [`FORMALIZATION_STATUS.md`](FORMALIZATION_STATUS.md). No Lean build,
Lean cache, or Lean verifier is required for the proof described here.
