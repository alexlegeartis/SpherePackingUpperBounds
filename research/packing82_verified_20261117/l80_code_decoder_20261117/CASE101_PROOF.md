# Case101: primary exact depth6 and81D extension (independent audit pending)

## Parent
Use the independently reconstructed Bachoc-Nebe neighbor in ../bachoc_nebe80_reconstruction_20261117. It is even unimodular of rank80, with explicit norm8 vector; minimum≥8 inherits Theorem4.7 of Bachoc-Nebe, Extremal lattices of minimum8 related to M22. The source's10-block E8, binary[10,5,4] code and neighboring hypotheses have been checked separately. No catalogue-isometry claim is needed.

## Coordinates and finite coverage
The parent basis has10-block E8 root coordinates M/2. The frozen half-coset is x=rM/4, where r is case101.json's binary vector. Its Gram norm rGr=2828 is divisible by4. Since G is even integral, every squared norm in x+L is an integer: (r-2v)G(r-2v)^T/4 = rGr/4-rGv^T+vGv^T. Thus excluding all distances≤5 proves depth≥6.

Every parent vector is uniquely described up to local lattice representatives by one of16 diagonal shifts z/2 plus an L' vector. Each L' vector has block coordinates s_i+2t_i. The s_i are residue representatives whose first4 polarization rows lie in C and whose last4 lie in Cperp. Both codes have dimension5 and are systematic on the first5 columns; the last5 states are exactly determined by the exported masks. All256^5 choices are covered by the DFS.

The neighboring condition is Σ(s_i+2t_i)∈2L1. Every codeword has even Hamming weight. Therefore the required XOR of the10 local E8/L1 carry labels of t_i is the4-bit vector whose jth entry is wt(the jth Cperp row)/2 mod2. The chosen diagonal representatives already generate the singular complement, so there is no additional κ offset.

With E twice the Euclidean E8 basis, X=8x in Euclidean coordinates, the local E8 target after division by2 is (X_i-4s_i E-2zE)/16. Its integer squared residual numerator q contributes q/128 to the actual80D squared distance. The global cutoff is640. Exact integer/half-integer parity rounding gives an E8 nearest point v; the code checks its distance≤1 for every local table. Any relevant local point y has distance²≤5/2 from that target, so ||y-v||<sqrt8. E8 norms are even, hence y-v belongs to the completely enumerated shell0,2,4,6 (9121 vectors). All local points relevant to cutoff640 are therefore covered. Minima are retained separately for all16 carry classes.

DFS lower bounds use unrestricted minima for unassigned systematic columns and relax unknown bits independently in the last5 columns. These can only reduce the possible completion cost. At a complete residue word, a16-state XOR dynamic program computes the exact minimum compatible carry cost. Both minima and sums are integers. The node cap produces NOT_CERTIFIED, never an absence claim. All16 glue choices must finish.

## Completed primary result
The rowwise-fill portable backend replay, case101_replay.json and progress log, completes all16 branches:5,877,535,660 DFS nodes and100 residue leaves, with no vector at squared distance≤5. This is primary evidence pending separately implemented completeness audit.32 zero-distance parent-lattice controls and64 independently checked positive rejection witnesses test the maps, but do not replace absence review.

The separate contact output supplies an exact parent vector v with (r-2v)G(r-2v)^T=24. build81.py checks its full80D membership directly. Hence the coset depth is exactly6 conditional on the complete exclusion proof.

##81D extension
Take a short representative h of this half-coset and generate <(L,0),(h,sqrt2)>. Odd layers have norm≥6+2=8; nonzero even layers have vertical norm≥8; the parent attains8. The Gram determinant is2. candidate81.json gives the full rational Gram and exact center-density square2^80. A rational alternating-series Machin enclosure proves the ordinary density floor2.87777416803990e-17. Compared with the pinned preceding repository bound this is+3.1054957484%; against the derived external MW80×A1 comparator it is+80.3894724750%. The former is an internal incremental comparison, not an external prior record.
