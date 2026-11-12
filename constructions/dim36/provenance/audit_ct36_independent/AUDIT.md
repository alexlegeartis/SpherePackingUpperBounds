# Independent 36D minimum audit

Independent reviewer: audit_36d_certificate.

Exact checks independently executed using SymPy:
- For t=1/5 and t=1/4, integer Gram B=5G4+20tD; reduced Q=UBU^T; det U=+/-1.
- Q is symmetric integral with even diagonal. Exact LDL identity Q=LDL^T has positive pivots.
- Every input interval encloses its exact rational LDL entry. Endpoints are integers less than 2^53 divided by 2^40, hence exactly representable binary64.
- Quarter candidate .npy agrees exactly with B=5G4+5D.
- Physical determinant at quarter: 10090298369529/4096. At fifth: 3340880266351419591675846359145849/10^24.
- e1 in original coordinates has integer Gram norm160 in both candidates, hence physical squared norm8.
- Density ratio to center density 2^18/3^10: quarter139968/117649, fifth2187000000000000/2140753641621841.

## Enumeration correctness
At recursion depth k, q is a lower bound for contributions already fixed. cl,ch enclose sum_{j>k} L[j,k]z[j]. Each multiplication and addition is expanded by nextafter; endpoints and integers are exactly represented. The radius upper bounds sqrt((159-q)/d_lower), with expansion after subtraction, division and sqrt. All true norm<=159 branches occur within ceil(-ch-radius) through floor(-cl+radius). The distance from zero to [x+cl,x+ch] is a lower bound for the absolute LDL linear term; downward-rounded squaring/multiplication and accumulated sum preserve the lower-bound invariant. Zero is not removed from recursion. The terminal exact __int128 form evaluation counts all leaves; one leaf must be zero. Integrality then gives minimum>=160, and the witness gives equality. Evenness is not essential when the searched threshold is159.

Unsafe large coordinates/radii abort the program, never prune silently. At allowed coordinate magnitudes <=10^6 and max |Qij|202, __int128 has ample headroom. Node counter has no budget/cap. Successful completion does not depend on a capped contact shell.

## Compiler/runtime
Independent rebuild: g++ -std=c++17 -O2 -fno-fast-math -ffp-contract=off -frounding-math.
GCC14.2.0, x86_64. Disassembly verifies scalar SSE binary64 multiply/divide and hardware sqrtsd. No unsafe reassociation/FMA. Correctly rounded IEEE arithmetic plus one nextafter expansion justifies enclosures on this platform. This is not a claim of portability to fast-math or non-IEEE targets.

Fifth independently completed108910425 nodes,1 leaf. Quarter independently completed118119521 nodes,1 leaf. Both reruns match producer counts.

Novelty/current best-known baseline must be audited separately; this audit certifies the supplied explicit lattice minimum and density, not exhaustiveness of current literature.
