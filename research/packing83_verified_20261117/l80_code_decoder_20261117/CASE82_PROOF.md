# Mixed-depth82D extension: primary complete, independent audit pending

Parent and case101 depth6 are already independently verified. Case108 uses the same frozen integer decoder at cutoff640 (distance5): all16glues complete,5,099,300,880nodes. Its target has rGr divisible by4; all half-coset distances are integers, so absence≤5 implies depth≥6.

Apply the exact block permutation recorded in case203.json's second_coset to case108. build82.py constructs the full80×80 integral automorphism A=M P M^-1, checks det±1 and A G A^T=G, and verifies the transported class. No assertion relies solely on the code's abstract automorphism group. The sum with case101 equals case203 modulo2L.

Case203's numerator norm is2 modulo4, so all half-coset distances are half-integers. The backend pair_depth55.cpp differs only by cutoff640→576; its local shell coverage remains valid at the smaller radius. All16glues exhaust2,489,483,280nodes with0residue leaves. Therefore the sum class has depth≥11/2, pending independent replay.

Let x,y represent the two depth6 classes. Attach vertical vectors a,b with Gram [[2,3/4],[3/4,2]]. Their form is Q(m,n)=2m²+(3/2)mn+2n². The least eigenvalue is5/4. The odd/even and even/odd parity classes have minimum2; odd/odd has minimum5/2, checked at(±1,±1), while all other odd pairs have squared coordinate length≥10. Nonzero even pairs have minimum8. Horizontal minima are respectively6,6,11/2, or0. Thus the full lattice has minimum8, attained in the parent. Its Schur determinant is55/16.

candidate82.json supplies the complete rational82Gram, exact automorphism, two half-coset numerators and comparison fractions. Exact center-density square is2^86/55; a rational Machin lower enclosure proves ordinary density≥8.56691236634380e-18. This is+14.0072495691% over the pinned previous repository result. The derived MW80×A2 example is an attributable external comparison, not asserted to be the strongest conceivable prior consequence.
