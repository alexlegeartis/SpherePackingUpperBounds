# Dimension-specific compact harmonic lattice bounds, 40–42D

Completed package for dimensions 40, 41 and 42. Each dimension has passed independent exact and interval audits plus a fresh full portable replay from this compressed package. See VERIFIED_RESULTS.json for the rigorous bounds and PORTABLE_REPLAY_REPORT.json for the release checks. No other dimension is certified by this package.

Install requirements, then run `python run_all.py` (two isolated dimension workers), or `python run_all.py --dimension 42 --jobs 1`. `--quick` checks exact algebra, independently implemented cap inequalities and the exact translated-polynomial tail, but is NOT a full geometry replay.

Full mode also replays all main geometry leaves at 384 bits, freshly encloses every scalar root, reruns all radial intervals and rechecks their exact budget and density normalization. Both tail arguments are retained: the original root-isolation proof and a supplemental exact proof that every coefficient of P(180+t) is negative. The scalar and main interval backends use the reviewed producer Arb/Jet arithmetic; independent checks cover exact linkage, cap algebra and complete trace traversal. See each dimension's theorem review for scope and hashes.

All files are hash-checked before use, lossless gzip is checked before and after decompression, and computation takes place only in a fresh temporary directory. Source files are not mutated. No result applies to unrestricted sphere packing.
