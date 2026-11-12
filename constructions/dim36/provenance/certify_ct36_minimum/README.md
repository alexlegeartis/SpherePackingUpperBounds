# Complete minimum enumeration for CT36 flex candidates

Independent minimum audit, 2026-11-12. No capped search is used as a proof.

## Inputs and normalization

`candidate_1.npy` in the sibling contact-walk directory is the integer Gram matrix
`Q0=5G4+4D=20 G(1/5)`. It is even integral. `candidate_quarter.npy` is
`5G4+5D=20 G(1/4)`, also even integral. To prove minimum squared length 8 in
physical normalization it suffices to exclude nonzero integer vectors of
Q-norm at most 158. Our exhaustive enclosure search uses bound 159, which
strictly contains every possible violation. A basis vector of Q-norm160
provides equality.

`reduce_bkz.py` applies BKZ to a rounded high-precision Cholesky embedding.
The approximation is used ONLY to find a change of basis: `prepare_exact.py`
checks det(U)=±1 and works thereafter with the exact integer Gram UQ0Uᵀ.
The selected transform is in `bkz36.npz`.

## Exact input and outward arithmetic

`prepare_exact.py` computes rational LDLᵀ and verifies its exact identity and
positive diagonal. Each diagonal and triangular coefficient is enclosed by
adjacent multiples of 2^-40, obtained by exact rational floor/ceiling.
The signed interval numerators fit in the exactly representable integer range
of binary64. The text input therefore encodes exact dyadic endpoints.

`interval_enum.cpp` enumerates from coordinate35 down to0. At each node,
its q is a LOWER bound on the already assigned LDL square contributions.
The interval center encloses the exact sum of triangular coefficients times
assigned integer coordinates. The radius is an UPPER bound computed using
q and the positive LOWER LDL diagonal. Every elementary interval endpoint
is rounded outward with nextafter; sqrt is rounded upward as well.
The new squared contribution uses the distance of the interval from zero,
thus a LOWER bound. Hence every vector of norm≤159 is visited, never pruned.
All leaves are checked by exact signed128-bit integer Gram arithmetic.
Runtime guards abort, rather than certify, for unsafe/nonfinite ranges.
No node cap, time cap, pruning heuristic, or approximate-radius reduction is
present. `COMPLETE 1` is emitted only after the entire recursion returns
without finding a violating vector.

Compiler command:

```
g++ -O3 -fno-fast-math -ffp-contract=off interval_enum.cpp -o /tmp/ct36_interval_enum
```

The implementation assumes ordinary correctly rounded IEEE754 binary64
operations and sqrt. With guard |coordinate|≤10^6, the exact norm accumulator
is safely within signed128-bit range for these small Gram coefficients.

## Results

At t=1/5: original exhaustive interval run completed108,910,425nodes,
exactly one leaf (zero). Independent exact input and algorithm audit by a
separate worker passed. Guarded repeat completed identically:108,910,425nodes, onezero leaf. Log is `guarded_interval.log`.
At t=1/4: guarded interval enumeration COMPLETE after118,119,521nodes,
exactly one leaf(zero), as recorded in `quarter_interval.log`. The exact
integer Gram has basis vectors of norm160, so physical minimum is exactly8.
Physical determinant is10090298369529/4096; the conditional comparison ratio
139968/117649 now has its minimum-distance premise certified. Novelty and
comparison against current literature remain separate checks.

Determinant, density formula, contact identities and baseline comparison are
separate obligations; see sibling `audit_ct36_flex`. This directory establishes
only the global minimum for each completed successful rigorous run.
