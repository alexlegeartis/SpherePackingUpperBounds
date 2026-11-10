# Independent audit: stronger dimension-four p4 candidate

PASS: fixed audited_snapshot.json yields density strictly below 0.634526023, subject to the published neighbor-count inputs. All five verification stages were rerun against an isolated snapshot, not a moving producer file. No defect found.

Snapshot: audited_snapshot.json.
SHA256: 3b8bba3d3bddccd178eec1527446a6462f7f49a9715a975721e6a48e7bd073c8.
Its parameters are R=333/200, J=1, p=4, K=16. Reconstructed exact matrices match the producer matrices.

## Reproduction

From this package directory run python verify.py. Equivalent component commands:
python certify_psd.py witness.json
python certify_radial.py witness.json
python certify_count_pair_arb.py witness_psd_exact.json
python certify_far_bernstein.py witness_psd_exact.json
python certify_cyclic_bernstein.py witness_psd_exact.json

Results stored alongside the snapshot:
- All 17 rational coefficient matrices positive definite by exact LDL. Smallest reported pivot: 9.900320832491374e-10.
- Fourier polynomial: zero positive roots, positive at zero.
- Physical tail: zero roots beyond conservative rational threshold below pi*R^2, negative at that threshold.
- Count-adjusted pair: 231 boxes, maximum depth 8, verified.
- Far: 135399 boxes, 33933 sign-passed and 33767 excluded, depth 36, verified.
- Cyclic: 127 boxes, 56 sign-passed and 8 excluded, depth 10, verified.

## Continuum and normalization review

The verifier source is unchanged from the previously audited 4D p3 candidate. The source/formula analysis and periodic-packing proof in ../INDEPENDENT_AUDIT.md apply with radial feature degree four. That analysis does not rely on degree three. In particular:

1. Exact rational matrices define a PSD kernel supported on {0} union the closed annulus [1,R]. The normalized angular basis U_k(t)/(k+1) is positive definite on S^3; higher angular features vanish at zero.
2. Gaussian-Laguerre radial functions are Schwartz, have Fourier eigenvalues (-1)^k, and both polynomial signs are global by exact Sturm sequences. The local Taylor verifier uses the same rational coefficients and matrices.
3. The cyclic Bernstein numerator is exactly (abc)^16 times the cyclic sum. It covers all sorted triples in [1,R]^3 including equality faces. Since R<2 all cube triples satisfy the triangle inequality. Symmetry covers the remaining permutations.
4. For exactly one long triangle side, only one kernel term remains; the far verifier covers it, including distance=R. Angular cosine is at most 1/2 in this regime, below the verifier's upper endpoint 1-1/(2R^2). Strict domain exclusions preserve boundaries. Two long sides make all terms zero.
5. The highest actual dyadic subdivision depth is 36, below binary64's exact endpoint precision limit. Outward rounding encloses each de Casteljau operation.
6. Right-hand pair subintervals include each count cutoff with the smaller allowance; this proves the half-open count condition. r=1 and R are included. Beyond R the physical tail handles the pair inequality.
7. Summing PSD local kernels gives nonnegative origin-plus-diagonal contribution because every unordered triangle contributes twice its nonpositive cyclic sum. Combining the count-adjusted pair inequality and Fourier Poisson bound yields density <= (pi^2/32)*(f(0)+K(0,0)+sum B_i lambda_i)/fhat(0). The support need not be continuous. Periodic replication of truncations and Schwartz decay extend to arbitrary packings.

The nonnegative count weights use bounds 24,25,26,27 at radii 1.008096,1.017445,1.025977,1.036950. The published input is https://arxiv.org/html/1206.2608 Section 6/Table 2: its rigorous dimension-four radii are each larger by 0.000001. Monotonicity therefore applies. This audit verifies applicability, not the original external count SDPs.

## Exact objective and updated baseline

Independent exact rational reconstruction:
- fhat(0) = 9999999999999148229243/10000000000000000000000.
- Unnormalized numerator = 205730968346242987111509968901/100000000000000000000000000000.
- Normalized numerator = 205730968346242987111509968901/99999999999991482292430000000.

256-bit Arb enclosure after multiplying by pi^2/32:
0.6345260220720713322566555674163184830382171356283089733216206914541048187119 +/- 4.17e-77.

Consequently density <0.634526023.

The updated comparison is arXiv:2512.04829v2, Table 1, https://arxiv.org/html/2512.04829, which reports 4D bound 0.63610277. Absolute decrease from that displayed bound is approximately 0.00157674792792866774334. Relative to the additive gap from that baseline to pi^2/16, the closed fraction is approximately 0.0818983686792133230, or 8.18983686792133230%. This is below the requested 20%. These numbers compare to the published displayed baseline, not its undisplayed exact value.

Scope: existing verifier source review, independent snapshot reruns, and independent exact objective reconstruction. This is not a clean-room implementation or formal verification of every calculation. The external neighbor-count certificates are published inputs. No unrestricted literature-priority claim follows from this audit alone.

The registry uses the outward-rounded bound 0.634526023 rather than the internal enclosure center. Its conservative gap decrease is 8.18983204813735% and upper-bound decrease is 0.24787614114619%. The portable witness.json contains exactly the same mathematical input fields as audited_snapshot.json, with optimizer provenance removed.
