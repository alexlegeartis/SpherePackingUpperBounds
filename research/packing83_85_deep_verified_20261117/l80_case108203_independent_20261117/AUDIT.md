# Independent mixed-coset 82D audit

PASS. Case 108 exhausts 5,099,300,880 nodes and 656 residue leaves; case 203 exhausts 2,489,483,280 nodes and zero residue leaves. All 16 glues complete. The separately authored integer decoder from the 81D audit is reused unchanged for 108 and with only CUT=576 instead of 640 for 203. Its shell enumeration, nearest-point coverage, carry DP and pruning proof remain valid at the smaller cutoff. No floating-point bounds occur. Node caps abort as INCOMPLETE.

The new map checker compares every non-target integer with the fully reconstructed 81D input, reconstructs both target arrays from the parent basis, and proves the entire transported-coordinate map is an integral unimodular Gram isometry. It checks that case 203 is the sum of case 101 and the transported case 108 modulo twice the lattice. Case 108 distances are integral and excluding <=5 yields >=6. Case 203 distances are half-integral and excluding <=4.5 yields >=5.5. Case 101 depth 6 is inherited from the previous completed independent audit.

The exact extension checker reconstructs all 82D Gram blocks and their Schur complement [[2,3/4],[3/4,2]]. Vertical parity minima are 2,2,5/2,8. Adding the respective horizontal class minima proves minimum 8; the parent attains it. Determinant is 55/16 and center-density square is 2^86/55. A rational Machin lower enclosure proves ordinary density >=8.56691236634380e-18.

Parent minimum 8 remains the Bachoc–Nebe Theorem 4.7 premise with previously checked constructive hypotheses. Exact arithmetic shares FLINT with the producer; the exclusion implementation is separately written. The external comparison is the pinned derived MW80 times rescaled A2 example, not an assertion of strongest prior literature.

Reproduce with PYTHONPATH=/tmp/packing-deps python check_maps.py and check_extension.py in the archived sibling layout. Compile search108.cpp and search203.cpp with g++ -O3 -std=c++17 and run against the matching case input.
