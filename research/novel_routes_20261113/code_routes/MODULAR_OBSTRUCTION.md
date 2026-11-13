# A universal obstruction to short-normal 71D sections

This is an exclusion result, not a new packing-density record. It was derived
jointly by the root and code-routes agents, with the Fourier-coefficient
calculation checked independently.

## Theorem

Let N be an even unimodular lattice of dimension 72 which has a nonzero vector
of squared norm at most 6. For every nonzero w in N with squared norm at most 10,
the section N intersect w-perp has a nonzero vector of squared norm at most 6.

Equivalently, if such a section has minimum at least 8, then its 72D parent
must itself have minimum at least 8.

The assertion applies to every even unimodular parent, including all prime
neighbors of the published extremal 72D lattice. It is not merely a statement
about the finitely many neighbors visited by a search.

## Proof

Put a=(w,w), and use the degree-two harmonic polynomial

    H_w(x) = (x,w)^2 - (x,x)*a/72.

The harmonic theta series

    Theta_(N,H_w)(q) = sum_(x in N) H_w(x) q^((x,x)/2)

is a level-one cusp form of weight 38. Write T_s for the sum of H_w over the
shell of squared norm s, so that its first coefficients are T_2, T_4, T_6.

The two-dimensional cusp-form space S_38 has basis

    F1 = Delta E4^5 E6 = q + 672 q^2 - 51084 q^3 + ...,
    F2 = Delta^2 E4^2 E6 = q^2 - 72 q^3 + ... .

Consequently every such theta series satisfies the exact relation

    T_6 + 72*T_4 + 2700*T_2 = 0.                         (1)

Suppose the section has no nonzero vector of squared norm at most 6. For
every x on any of the three short shells, integrality gives (x,w) a nonzero
integer. Hence

    H_w(x) >= 1 - (x,x)*a/72 >= 1 - 60/72 = 1/6 > 0.

All three T_s are therefore nonnegative, and at least one is positive since
the parent has a short vector. This contradicts (1), whose three coefficients
are positive. This proves the theorem.

## Consequence for the proposed neighbor construction

For primitive v of squared norm 24 in an extremal even unimodular Gamma72,
consider the even 2-neighbor

    H + Z*(v/2),   H={x in Gamma72 : (x,v) is even}.

It contains the squared-norm 6 vector w=v/2. The tempting construction was
to arrange that its section by w-perp had minimum 8 and determinant 6. Its
center density would be 2^(71/2)/sqrt(6), about 0.885% above the verified
norm 10 two-translate antipode packing. The theorem rules this out for every
v, without any enumeration.

In fiber language, there must exist x in Gamma72 with

    (x,v)=12,   (x,x) in {8,10,12}.

Indeed a short section vector in the odd neighbor coset has the form
x-v/2, has squared norm x^2-6, and gives one of those three possibilities.
The old shared sublattice has minimum 8, so the short section vector cannot
come from it.

The same theorem rules out the analogous minimum 8 sections with normal
norm 2 or4 in any even unimodular 72D parent. For normal norm 8 or 10 it says
that switching to a parent with minimum below8 cannot produce a minimum 8
section. It does not exclude rational metrics, odd unimodular parents,
periodic packings, or other section dimensions.

## Replay and standard mathematical inputs

Run:

    python research/novel_routes_20261113/code_routes/verify_modular_obstruction.py

The checker uses rational/integer arithmetic to compute the Fourier
coefficients from E4, E6 and Delta=(E4^3-E6^2)/1728; it checks the basis
dimension, relation (1), and the positive margin. It does not mechanize the
classical harmonic-theta modularity theorem.

Primary sources for the standard inputs:

- Scott D. Kominers, *Configurations of Extremal Even Unimodular Lattices*,
  International Journal of Number Theory 5 (2009), 457–464, section 2.2;
  <https://arxiv.org/abs/0706.3082>.
- William Stein, *Modular Forms: A Computational Approach*, chapter on
  level-one modular forms;
  <https://wstein.org/books/modform/modform/level_one.html>.

The existence and minimum of the extremal 72D parent are inherited from
Nebe's construction; they are not needed for the general theorem itself.
