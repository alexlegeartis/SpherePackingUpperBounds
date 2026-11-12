# A verified 36-dimensional lattice packing

Public research release, 12 November 2026. Computer-assisted verification package; not a peer-reviewed result.

## Result

The explicit rational Gram matrix `data/Q_original.txt` divided by **20** defines a positive-definite lattice of minimum squared length **8** and Gram determinant **10090298369529/4096**. Its center density is **16777216/3176523 = 5.281629001...**. Its ordinary sphere-packing density is this number multiplied by **pi^18/18!**, approximately **7.33035e-7**. A Gram determinant is the square of the covolume, not the covolume itself.

The KP36 reference has minimum squared length 8 and Gram determinant 3^20, hence center density 262144/59049. The exact ratio is **139968/117649**, or **18.970836981190% higher**. This is a packing-density lower bound, not a kissing-number result and not a proof of optimality.

## Idea and exact construction

The Coxeter–Todd/Turyn 36D seed admits a three-block metric deformation. Optimizing only its three scale parameters reaches a vertex with minimum 8 but still fails to beat KP36. Its known contact vectors leave an off-diagonal symmetric deformation direction D: their quadratic forms v^T D v vanish. Moving in this direction decreases the determinant without immediately shortening the contacts. Numerical null-space discovery is used only to suggest D, which is saved as an exact integer matrix.

Take G=G4/4 from `data/G4.txt` and D from `data/D.txt`. The candidate is G+D/4, equivalently Q/20 where Q=5G4+5D. Exact algebra gives

    det(G+tD) = 3^18 (5/4)^12 (1-16t^2/15)^12.

The selected parameter is t=1/4. A contact-preserving direction is NOT alone a minimum proof: previously unknown vectors could become short. The complete enumeration below closes that gap for this particular rational candidate.

The saved inputs fully specify the construction without optimization, external downloads or randomized basis reduction. `generate_and_check.py` reconstructs it from G4 and D, verifies the supplied unimodular transformation U, positive rational LDL factorization, determinant, density comparison and exact norm-160 witness e1. It generates the complete enumeration inputs and matches the historical audited inputs exactly. It does not itself prove global minimum.

## Minimum proof and reproduction

Requirements: Python 3 with SymPy (tested 1.14.0), GCC C++17 on x86_64 with correctly rounded IEEE754 binary64 arithmetic and sqrt. No fast-math. NumPy is needed only for the historical discovery scripts, not the standalone verification.

Run `python3 verify_all.py` for integrity checks, exact reconstruction, compilation and the complete minimum replay. To run the stages separately, from this directory:

```sh
python3 generate_and_check.py
# Independent audit's build configuration:
g++ -std=c++17 -O2 -fno-fast-math -ffp-contract=off -frounding-math interval_enum.cpp -o verify_minimum
./verify_minimum ldl.txt Q_reduced.txt > minimum.log 2> progress.log
cat minimum.log
```

Expected output:

    COMPLETE 1 nodes 118119521 leaves 1

The integral Gram Q is even; a nonzero physical norm below 8 would have Q-norm at most 158. Searching through threshold159 includes every possible violation. Exact rational LDL coefficients are enclosed by exact dyadic intervals with denominator2^40. The recursive enumerator bounds every center and radius outward using `nextafter`; accumulated squared contributions are bounded below. It therefore cannot prune any vector of Q-norm at most159 under the stated floating arithmetic assumptions. All leaves are evaluated with exact signed128-bit integer arithmetic. Only the zero leaf occurs. Unsafe ranges abort rather than certify. A basis vector with Q-norm160 gives equality.

**This proof is computer-assisted and depends on correctly rounded binary64 operations.** It is not a wholly symbolic enumeration or a claim of portability under arbitrary compiler options. An independent reviewer audited exact input preparation and the recursion, rebuilt with the above independent flags and reran the entire enumeration, obtaining the same counts. This is an independent rerun and source audit of the same implementation, not a separately designed second enumeration algorithm. The original guarded build used `-O3 -fno-fast-math -ffp-contract=off`; see saved logs and audit.

## Contents and provenance

- `data/`: plain-text exact G4, D, U and candidate Q.
- `generate_and_check.py`, `interval_enum.cpp`: self-contained construction and verification entrypoints.
- `exact_checks.log`, `environment.txt`: package preparation checks and environment.
- `provenance/certify_ct36_minimum/`: original verifier, interval inputs, basis data and complete runs, including the earlier t=1/5 result.
- `provenance/audit_ct36_independent/`: independent audit, hashes and full rerun outputs.
- `provenance/audit_ct36_flex/`: exact determinant/deformation checks and exploration.
- Other provenance directories preserve the seed, generation scripts, transformations and partial contact shells.

Historical READMEs deliberately retain their original conditional or 'unverified' labels. They describe the state before the final complete minimum verification and are not the current status of t=1/4. Historical exploratory scripts may contain original paths; they are provenance, not the portable verification entrypoint. Their bounded floating searches and partial shells are not used as minimum proofs. No unfinished stronger candidates are claimed here.

## Sources and comparison scope

1. Nebe–Sloane lattice catalogue, CTturyn2 seed: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/CTturyn2.html . Original source HTML is retained under provenance/ct36_metric_screen_20261112/source.html.
2. Nebe–Sloane lattice catalogue, KP36 baseline: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/KP36.html . Minimum8, determinant3^20.
3. Jumagulov, July2026, arXiv:2607.11319v2: https://arxiv.org/abs/2607.11319 . The project's literature check reports KP36 as best known36D packing, not just best lattice.
4. Dorofeev–Sun–Wang, August2026, arXiv:2607.20359v4: https://arxiv.org/html/2607.20359v4 . The project's check found records38–43, not36.
5. Nebe's CT polarization context, section5.2: https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/Ext72.pdf . This is context, not a proof of this deformation.

The separate literature check found no later36D improvement. That search is not an exhaustive novelty guarantee. Direct live access to Cohn's table failed; this package does not claim that table was checked successfully. The exact mathematical comparison with the specified KP36 baseline does not depend on the completeness of the literature search.

SHA256SUMS inventories all package files except itself. No compiled executable is needed or trusted: rebuild from source.
