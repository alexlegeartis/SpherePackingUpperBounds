# Independent 83D exclusion and density audit

PASS. The independently authored integer decoder from the completed 81D audit was reused with only CUT changed to 384 (case 301) and 512 (case 302). Complete fresh traversals exhausted all 16 glues: 7,416,859 nodes and 150,439,952 nodes respectively, both zero residue leaves. The second node count differs from the primary because the row/tie ordering differs. No floating pruning occurs; computation limits throw INCOMPLETE.

check_maps.py independently reconstructs each target from the parent basis, compares every common decoder array with the previously reconstructed input, and verifies numerator norm modulo 4 is zero. Thus all distances are integral: excluding <=3 and <=4 proves depths >=4 and >=5. The previous local completeness proof at cutoff 640 remains valid for these smaller cutoffs.

The separate l80_83_geometry_independent_20261117 audit reconstructs all five integral unimodular isometries, all class relations, eight parity classes and their global vertical minima. I reviewed its exact spectral lower bound and finite cube argument and freshly replayed it. Parent minimum and cases 101,108,203 are inherited from completed independent audits, not new unverified premises.

check_density.py compares the full primary Gram to the independently reconstructed Gram, verifies determinant 23/4, center-density square 2^85/23 and a rational Machin lower bound giving density >=2.56962386473426e-18. External comparison is the pinned derived MW80 x A3 construction, with A3 rescaled from minimum2 to14 and determinant4*7^3. This is not a claim of strongest prior literature.

The arithmetic backend FLINT is shared; the exclusion implementation and geometry implementation were written separately from the producer. No new independent claim is made for the underlying Bachoc–Nebe source theorem.
