# Independent review of the Ne64 coset interval enumerator

Reviewed by the code-routes research agent before the complete four-shard
run. This is a mathematical/code audit; the actual completion evidence
and exact65-dimensional Gram matrix belong to `../ne64_hole_certificate/`.

Files reviewed: `prepare.py` and `interval_enum.cpp` in that directory.

## Exact reduction and congruence

The setup verifies `Q=T G T^t` with integral unimodular T and an exact
inverse, and checks the candidate residue relation `cT=a mod3`. It also
checks `aG=0 mod3`, making a/3 a dual vector, and `aGa^t=6 mod18`.
All representatives of a/3+L therefore have norm congruent to2/3 modulo2.
Their integer numerator norms are6 mod18, so excluding numerator norms
at most60 proves minimum at least78/9=26/3. The setup's exact LDL
factorization also proves positive definiteness.

## Why the pruning is exhaustive

For row coordinates z and `Q=L D L^t`, the norm is

```
sum_i D_i (z_i + sum_{j>i} z_j L_{j,i})^2.
```

The preparation file writes precisely these column coefficients, and the
C++ enumerator uses them in that orientation.

Each coefficient is enclosed between exact integer multiples of2^-40.
The upper-coordinate products and sums are accumulated in signed128-bit
integers before conversion to binary64 and one outward `nextafter` step.
The enforced limits give numerator magnitude below2^79, so this integer
accumulation cannot overflow. Scaling by2^-40 is exact in the relevant
finite range.

The recursion carries a lower bound on the assigned norm. Subtracting it
from60 with upward rounding and dividing by a lower diagonal bound gives
an upper bound for the remaining squared radius. Outward-rounded square
root and integer endpoints consequently contain every possible next
coordinate. The interval's minimum absolute value gives a lower bound
for that coordinate's cost. The subsequent products and sum round
downward. Thus every pruning test discards only a branch whose true norm
exceeds60.

The coordinate cap and Gram-entry cap also bound each exact leaf norm
well inside signed128-bit range. A cap violation aborts rather than
silently discarding a branch. Compilation must retain `-fno-fast-math`
and `-ffp-contract=off`, as documented in the source.

## Shards

The four shards enumerate the same upper-coordinate prefixes in the same
order, assigning each prefix by its index modulo4. Lower-coordinate
values left by a completed branch cannot affect a later prefix, because
every centre depends only on already assigned higher coordinates.
An accepted completion therefore requires all four shards to report
COMPLETE with found=false, identical total-prefix counts, assigned-prefix
counts summing to that total, and input hashes matching the exact setup.

The coordinating agent subsequently reported all four completions:
4,880,529,145 total visited nodes and42,255,279 common prefixes, with no
leaves or witnesses. See the actual certificate summary for the saved
records; this note does not replace those records.
