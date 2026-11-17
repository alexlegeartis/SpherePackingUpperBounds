# Recursive binary-offset flags and the 95D improvement

This extends the 94D parallelogram proof. The parent and covering constant are unchanged.

Choose horizontal generators a_j successively, for j=1,2,3. At step j exclude the 3^(j−1) congruent torus balls centered at all sums of earlier generators with coefficients in {−1,0,1}, each of squared radius d_(3^(j−1)). The strict covering inequality leaves a point a_j outside these balls. This guarantees that every signed sum whose highest nonzero coefficient is at position j has squared torus norm exceeding that depth.

Assign the eight subset sums of a_1,a_2,a_3 to the eight coordinate parity colors of a rank-three vertical lattice. If the highest differing color bit is j, the horizontal difference has the required depth d_(3^(j−1)).

Write h_j=16−d_(3^(j−1)), so 4<h_1<h_2<h_3. Use the vertical Gram

    [[h1, −h1/2, 0],
     [−h1/2, h2, −h2/2],
     [0, −h2/2, h3]].

Its quadratic form is

    h1 A_3(z1,z2,z3)/2
    +(h2−h1) A_2(z2,z3)/2
    +(h3−h2) z3²,

where A_k denotes the standard root-lattice quadratic form. Each A_k/2 is an integer, is nonnegative, and is at least1 on a nonzero integer vector. Thus if the highest nonzero coordinate is at position j, the vertical norm is at least h_j. The highest nonzero coordinate is at least the highest odd coordinate, proving every nonzero parity class has the required lower bound. Nonzero even vectors have norm at least4h1>16. The combined periodic packing therefore has squared minimum16.

The period is L orthogonally summed with2M and has eight translates. Center density squared is

    (2^184/47^14) · 64/det(M).

The determinant is h1 h2 h3 − h1 h2²/4 − h1² h3/4. The exact arithmetic check flag95.py verifies the occupied-density floor

    2.06770138493108 × 10^−21.

For odd dimension95 the squared unit-ball volume is exactly 2^96 π^94/(95!!)^2. The saved rational Machin lower bound gives a rigorous lower density, independently of the decimal mpmath display. This improves the prior displayed bound by approximately1.20187%.

The same weighted-chain proof works in arbitrary extension rank, but the four-dimensional A4 chain does not beat the existing D4-based 96D construction. No 96D improvement is claimed.
