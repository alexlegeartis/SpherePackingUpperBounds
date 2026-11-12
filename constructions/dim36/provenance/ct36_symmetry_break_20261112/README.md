# New high-dimensional construction screens, 2026-11-12

**No new packing record.** This folder contains a fresh finite-code rate trade, rigorous cross-dimensional pruning, and a cheap metric/fiber screen. Benchmarks come from sibling density_lower_baseline/records_1_48.json, general (not lattice-only) records.

## Highest-priority new idea: pay ternary redundancy to buy quaternary rate

The exact KP mixed hexagonal quotient from sibling kp36_code has local squared costs0 at(0,0),4 at(0,±1),3 at(nonzero,0),1 at(nonzero,±1), and period norm12/covolume6sqrt3. Use an outer quaternary code C of length18, minimum6, dimension11 and a ternary code ker(H), where H has rank3 and no zero column.

A sufficient condition for minimum8 is: for every support S of a weight6 or7 word of C, no choice of signs ε_i∈{±1} satisfies sum_{i∈S}ε_i H_i=0. Proof: any zero trit on S costs two more than its cheapest value; any nonzero trit outside S costs four more. Thus below8 is possible only with exactly S as the nonzero trit support, which this condition excludes. Same-outer differences have ternary weight≥2 because H has no zero columns. Words of weight≥8 are automatically safe.

If feasible, total point count relative to KP36 is4²/3²=16/9: **+77.7778% in36D**. This is a CONDITIONAL target, not a construction yet. The rank2 version with a dimension10 outer code would give4/3. Full fiber parity-phase changes previously tested are different: this proposal genuinely reduces ternary fiber cardinality to eliminate short support patterns.

`outer18.py` reconstructs Grassl's [18,11,6] by shortening the provided [20,13,6] generator at positions19,20. Complete exact enumeration of4,194,304 words gives1467 weight6 and3804 weight7 words, with1757 distinct supports. Files outer18_G11.npy and outer18_bad_supports.npy are the actual finite feasibility instance. No H has been found or certified here. Audit/feasibility assignment was handed to the parallel audit worker.

`check_patterns.py` rigorously classifies rank2 nonzero projective columns: a length6 support avoids a full-support signed dependency iff its color-count multiset is(5,1,0,0); for length7 it is(6,1,0,0). Rank3 remains the more promising target. The familiar four projective columns of F3² have no full-support signed dependency at length4, illustrating that signed-sum avoidance is not vacuous, but not solving the length6/7 system.

Source: https://www.codetables.de/BKLC/BKLC.php?q=4&n=18&k=11 ; its cited source is Edel–Bierbrauer, Inverting Construction Y1 (1997). Local packing construction source and exact proof are in sibling kp36_code/README.md.

## General nonlinear distance-eight outer-code sweep

For length n and arbitrary quaternary code size M, minimum distance8, with a full ternary parity fiber, the2n-dimensional center density is M/3^(n/2+1). `code_sweep.py` generates rationally verified Delsarte dual polynomials. It proves this entire family cannot improve current records in26,28,42,44,46,48D. Remaining targets:

|dimension|outer length|minimum M for improvement|rigorous LP ceiling, rounded down|
|--|--|--|--|
|30|15|11364|15536|
|32|16|50503|53798|
|34|17|75754|179175|
|36|18|262145|661823|
|38|19|969645|2526885|
|40|20|7272337|8612913|

This prioritizes new nonlinear structures in34–38D; it does not assert existence below the ceiling. Shortening the known18D code gives65536 words in length17, so34D needs15.59% more rather than a one-word augmentation.

`puncture_spectrum.py`: integer Walsh spectra plus Hoffman bound prove **all18 single-coordinate puncturings** of the known length18 code have maximum distance-eight subset≤65536, below the34D target. Twelve sampled coordinate-pair puncturings have bound1736704/37≈46937.946, below the32D target; this is not all153 pairs. Twelve sampled triple puncturings leave sufficient spectral headroom, so `quotient30.py` tests three triples and every15 possible extra shortening coordinate. All45 exact256-node Cayley clique searches give4 cosets×1024=4096 words, below11364. No claim about arbitrary subsets or other codes follows from this quotient restriction.

## CT36 symmetry breaking and torus addition screen

Starting at the previously derived best three-block vertex(a,b,c)=(1/2,1,5/2), the Gram determinant is3^18*(5/4)^12. Explicit norm8 vectors bound its density by12288/15625=78.6432% of the36D record. Full metric symmetry breaking, unlike the already excluded three-scalar family, was not globally excluded.

A floating shell enumeration stopped at100million nodes and yielded35864 explicit norm8 vectors after LLL. Every listed vector has exact norm8 under the rational Gram. **Enumeration was incomplete**, so it is not a minimum certificate. A4000-vector floating eutaxy feasibility LP timed out after25seconds, which proves neither eutaxy nor non-eutaxy.

A second-copy translate could in principle beat the record if its squared cross-distance exceeds8*(15625/24576)^(1/18), approximately7.80. We sampled900 rational shifts (300 each denominator2,3,4). The cheap coordinate decoder left16 apparently viable candidates; descent using actual norm8 vectors rejected all16. `verify.py` recomputes all900 integer nearest-point witnesses and verifies their two-copy density upper bounds are≤the record. This is a sampled rejection, not a global covering-radius result. It demonstrates why a high approximate hole estimate must not be called a packing.

## Reproduction and verification

Run `python research/ct36_symmetry_break_20261112/verify.py` to independently check the12 rational dual certificates and900 rejection witnesses. Other named scripts regenerate the finite searches. `enum.cpp` is exploratory floating enumeration and has an explicit100million-node cap. It must not be used as an exact global minimum verifier.
