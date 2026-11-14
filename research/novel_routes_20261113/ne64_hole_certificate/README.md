# A certified Ne64 dual coset of minimum 26/3

All four exhaustive interval shards completed. This certifies an explicit
dual coset of `Ne64` with minimum **26/3**, and hence the explicit 65D lattice
of minimum 12 and determinant `10*3^31` in `../higher_parent/`.

## Exact reduction to a finite computation

Let `G` be the published 64-dimensional Gram, and let `Q=T G T^t` be the
saved integral Gram after an exactly verified unimodular basis change.
The selected dual coset has row coordinates `h=a/3` in the original basis.
The preparation script checks `aG=0 (mod 3)`, `a != 0 (mod 3)`, and
`aGa^t=6 (mod 18)`. Every numerator in that coset therefore has squared
norm congruent to 6 modulo 18.

Writing `c=aT^-1 (mod 3)`, the required computation is

\[
\{z\in\mathbb Z^{64}:z\equiv c\pmod3,\ zQz^t\le60\}=\varnothing.
\]

The next allowable numerator norm is 78, so emptiness proves that the
coset minimum is at least `78/9=26/3`. The exact norm-78 witness saved in
`../higher_parent/ne65_candidate.npz` proves equality.

## Why the enumeration encloses every relevant vector

`prepare.py` computes a rational LDL decomposition
`Q=L diag(d) L^t`, verifies it exactly, and encloses every coefficient in
dyadic intervals with denominator `2^40`. The binary64 endpoints are
exactly representable.

The recursion assigns coordinates from 63 down to 0. At each node it:

1. Accumulates lower and upper center numerators exactly in signed
   128-bit integers, then converts them outward to binary64 intervals.
2. Maintains a lower bound on the already assigned squared norm.
3. Computes an upper bound for the remaining coordinate radius, using
   the lower diagonal endpoint and outward-rounded arithmetic.
4. Visits every integer in the resulting interval with the prescribed
   residue modulo 3.
5. Updates the partial lower bound using the distance of the center
   interval from zero. Any leaf is checked using the exact integer Gram.

Every floating operation involved in exclusion is enclosed outward with
`nextafter`. The build disables unsafe fast-math and contraction. Input
and coordinate caps keep all signed integer operations within range;
an unsafe bound causes an abort, not a successful certificate.

The tree is partitioned after assigning the top 16 coordinates. Each
shard numbers the same prefixes in the same deterministic order and
handles precisely one residue class of prefix indices modulo four.

## Recorded completion

| Shard | Nodes | Assigned prefixes | Status |
|---|---:|---:|---|
| 0 | 1,218,700,383 | 10,563,820 | COMPLETE |
| 1 | 1,216,139,957 | 10,563,820 | COMPLETE |
| 2 | 1,219,964,465 | 10,563,820 | COMPLETE |
| 3 | 1,225,724,340 | 10,563,819 | COMPLETE |

All shards report the same **42,255,279** shared prefixes; assigned counts
sum to that total. The combined **4,880,529,145** nodes include repeated
upper-tree work. No shard reached a leaf or found a vector of norm at
most 60. The observed run took approximately three minutes of wall time
on four cores after exact preparation. Runtime on another machine may
differ.

`summary.json` contains the complete results, input hashes, source hashes,
and arithmetic assumptions. `run_manifest.json` records the compiler,
flags, input hashes and binary hash. `summarize.py` checks the complete
prefix partition and all four completion records.

## Reproduce

From the repository root, with Python, NumPy, SymPy, and g++ available:

```bash
OPENBLAS_NUM_THREADS=1 python research/novel_routes_20261113/ne64_hole_certificate/run.py
```

This reruns exact preparation, compiles the interval enumerator, executes
all four shards, aggregates them, and runs the exact 65D lattice verifier.
The complete traversal is the certificate computation; checking a stored
JSON status alone does not reproduce the exhaustive proof.

The parent minimum 12 is inherited from Nebe's published theorem,
as recorded in the [Ne64 catalogue entry](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CQ64b.html).
The exact catalogue Gram is retained and independently matched. This
certificate makes no assertion about global sphere-packing priority.
