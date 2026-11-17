# Three overlapping deep planes in the L80 parity quotient

This supersedes the preliminary candidate.json (det45/8); use candidate_alt.json.

The pinned rank80 minimum8 unimodular lattice has Gram G and twice-root-coordinate basis M in ../bachoc_nebe80_reconstruction_20261117/lattice80.npz.
Choose half-coset numerators r1,r2,r3 from candidate_alt.json. For six masks1..6, `known` records exact common E8 isometry and block permutation transferring from source101/108, whose independent depth6 audits are inherited. All six half-cosets thus have squared distance at least6 from the parent lattice.

For mask7, alt7.txt is the same frozen decoder input as case101.txt except target X=(rM)E, r=r1+r2+r3 mod2. The exhaustive code-derived E8 local decoder at cutoff384/128=3 returns no vector for all16glue classes (849936nodes,0leaves). Since rGr is divisible4, all squared norms in r/2+L are integers. Therefore its distance is at least4. This fresh decoder claim requires independent replay before publication.

The extension has generators (basis(L),0) and (ri/2,vi) for i=1..3, where vertical Gram V has diagonal2 and offdiagonal1. For an integer triple z, its horizontal coset is determined by z mod2. The minimum vertical squared norm is2 for every nonzero parity except111, where it is4. Zero parity has nonzero minimum8. Thus every nonzero extension vector has squared norm at least8. The matrix V=I+J has determinant4 and eigenvalues1,1,4; finite enumeration [-2,2]^3 plus ||z||²>8 tail suffices to check the parity minima.

The full lattice Gram is [[G,G R^t],[R G,R G R^t+V]], R rows ri/2, determinant detG detV=4. Its center density squared is2^83/4=2^81. The prior83 extension determinant23/4 is improved by ratio sqrt(23/16)=1.1989578808... . No claim of global priority without external comparison.
