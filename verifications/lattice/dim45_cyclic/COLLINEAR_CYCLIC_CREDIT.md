# Exact lattice-only collinear residual extraction

Use the hypotheses and notation of CYCLIC_CLASS_BRIDGE.md. The explicit
published kernel depends only on the two radii and their mutual angle, so
K(-x,-y)=K(x,y). Define, for v!=0,
 R(v)=K(v,-v)+2K(v,2v).
This is precisely the cyclic sum for the three-point configuration{-v,0,v},
so R(v)<=0. It is not a standalone pairwise sign assumption.

Inside the ordered off-diagonal sum O, take the subset of pairs (x,y) with
 y=-x, y=2x, or x=2y. The three cases are disjoint for x,y nonzero. This subset
is a union of complete three-anchor cyclic orbits (including both orientations).
Its sum is exactly sum_{v!=0}[K(v,-v)+2K(v,2v)]=sum R(v).
The complementary subset is also a union of complete cyclic orbits, hence
its sum is<=0. Absolute Gaussian convergence justifies all reindexing.
Thus O<=sum R(v). PSD gives
 C+sum D(v)+O>=0, hence sum[D(v)+R(v)]>=-C.

Therefore replace D by D_eff=D+R in the hybrid spatial majorant, leaving the
objective penalty +C unchanged. For t=|v|^2 and the published kernel,
 R(t)=C/B^2 [h(t)^2 Q(-1)+2 B h(t)h(4t)],
 D_eff(t)=C/B [(1+Q(-1)/B)h(t)^2-2h(t)+2h(t)h(4t)].
This formula is per vector. A charge written per opposite pair{+v,-v} has an
additional factor2; class-budget vector sums already supply that factor.
Do not add it twice.

For frozen38D Q, Q(-1)/B is approximately -0.13098055864982353. At t=1,
h(1)=1 and h(4)=4^96 exp(-360) is tiny but positive and must be retained or
bounded rigorously. Thus the additional negative correction is roughly13.1%
of the original kernel correction at t=1. This is a local correction size,
NOT a demonstrated density improvement. New radial envelopes and the final
profile certificate remain necessary.
