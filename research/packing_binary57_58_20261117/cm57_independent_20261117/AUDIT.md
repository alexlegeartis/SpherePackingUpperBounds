# Independent MW56 coset/57D/58D audit: PASS

Completed: all frozen hashes; independent rational LDL identity and every serialized binary64 enclosure; exact unimodular reduced basis; evenness and modulo6 congruence; contact numerator norm60 and original coordinate transport; full57D/58D block Grams, determinant and center density; fresh140-section curve/height/positive-definiteness/generation replay. Separate200/201-term rational Machin lower bounds verify floors57D1.53307111021240e-11 and58D5.02419365648403e-12 conditional on the completed absence proof.

The parent is minimum10, not14. The characteristic29 curve has30 simple finite zeros of1+t^30, all irreducible typeII; infinity is smooth after x=t^10X,y=t^15Y. Therefore chi=5. Schütt–Shioda's explicit height formula (https://arxiv.org/html/0907.0298 , Section11.8) gives h(P)=10+2(P.O)>=10 for nonzero sections; integral sections attain10. Section11.17 explains the degree/integrality criterion. The140-section verifier recomputes all curve equations and pairwise heights and proves that their coordinates generate the saved Gram lattice. These are explicit arithmetic checks combined with that external geometric height theorem, not a brute-force parent shortest-vector enumeration.

The producer's binary64 verifier was independently reviewed: lower partial quadratic contribution, enclosing centers, outward remaining radius and all enclosing integers, abort on overflow/nonfinite/caps. Its fixed2^-10 center uncertainty is safe under its input and coordinate restrictions: coefficient uncertainty<2^-26/term, arithmetic terms<2^-28/term, accumulated additions<2^-14, total<2^-13. Leaf128-bit evaluation cannot overflow under the checked bounds. Exact input enclosures must be checked separately, as done by check_setup.py.

A separately written exhaustive solver independent_search.cpp completed successfully. It uses long double (at least64-bit significand), ascending enumeration instead of center-out ordering, and a tighter independently derived fixed center enclosure2^-16. The coefficient uncertainty contributes<2^-19 total; long-double arithmetic contributes<2^-24, so this enclosure remains conservative. It uses nextafterl-style outward bounds, exact128-bit leaf values,6b node cap and no pruning heuristic. Input is the same hash-checked binary enclosure file. Basic IEEE operations and correctly rounded sqrt remain platform assumptions. It is not an integer-only or formal proof.

Build: g++ -O3 -std=c++17 -ffp-contract=off -fno-fast-math -frounding-math independent_search.cpp -o independent_search
Run: ./independent_search input/random56_37_bits.txt

The completed independent log says EXHAUSTIVE_NO_VECTOR, with3647514565 nodes,0 candidate leaves and969.20 seconds. It did not hit the6b cap. This is separate from the producer's3791478090-node binary64 report. A control build differing only in the admissible input dimensions passed40 exact brute-force tests,23 with witnesses and17 empty. Tests support implementation review but do not replace the mathematical enclosure argument.


## Construction closure

Every numerator norm in this class is divisible by6: (r-3z)G(r-3z)^T=rGr^T-6zGr^T+9zGz^T, with even G and rGr^T divisible by6. Excluding all numerators<=54 therefore proves norm at least60/9. The independently verified contact numerator has norm60, so the exact coset minimum is20/3. The opposite coset has the same minimum by negation, and3x belongs to the parent lattice.

For57D add(x,sqrt(10/3)). Nonzero layers whose index is not divisible by3 have horizontal norm>=20/3 and vertical norm>=10/3. Nonzero indices divisible by3 have vertical norm>=30; the zero layer has parent minimum10. Thus minimum10 is exact and determinant(10/3)*29^8.

For58D add(x,a),(-x,b), with vertical Gram[[10/3,5/3],[5/3,10/3]]. A layer(u,v) withu-v not divisible by3 has horizontal norm>=20/3 and vertical norm>=10/3. When u-v=0mod3 and(u,v) is nonzero, u²+uv+v² is a positive multiple of3, giving vertical norm>=10. The zero layer has parent minimum10. The verified block Gram therefore has minimum10 and determinant(25/3)*29^8.

The exact center-density squares are(5/2)^n/determinant. A separately coded200/201-term rational Machin lower pi bound proves ordinary density floors57D1.53307111021240e-11 and58D5.02419365648403e-12 by squared inequalities, including odd-dimensional normalization. No decimal approximation is used to prove these floors.

All mathematical and arithmetic checks passed. This is a new complete independent implementation/replay, not formal verification or an integer-only enumeration. Runtime IEEE rounding and correctly rounded elementary square root assumptions are explicit. No exhaustive external-priority review or permission to publish is established by this audit.
