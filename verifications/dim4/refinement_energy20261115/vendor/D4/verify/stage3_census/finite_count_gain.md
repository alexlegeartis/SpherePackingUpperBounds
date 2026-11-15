# Recovering finite-cardinality overcharge in the saved spherical proofs

This is a written algebraic consequence of the three spherical certificates
checked in stage 2. It changes no coefficient or Euclidean sign. `run.py` in
this folder checks its inputs and recomputes its arithmetic; Section 5 of the
paper states the same argument.

## Exact multiplicities

For a code C of N unit vectors, write H_g(C)=sum_(x,y in C) g(x.y).
The original, unrepaired anchor energy at u is S_u, with origin kappa,
origin coupling v and kernel K. Its positivity is supplied by the saved
integer-Gram representation, unrestricted for the bounded-N26 certificate and
restricted to the declared N for either fixed-N certificate. Put

    P(a) = g(a)+2v(a)+K(a,a,1),
    C(a,b,t) = K(a,b,t)+K(a,t,b)+K(b,t,a),
    B_raw = kappa+g(1).

The exact identity is

    sum_u S_u = N*B_raw-H_g(C)
                + sum_(x!=y) P(x.y)
                + 2 sum_{unordered triples} C(a,b,t).

There are N(N-1) ordered pairs and binomial(N,3) unordered triples. The saved
checkers certify P<=eta and C<=3epsilon on their entire closed domains; this
is precisely what their subtraction of eta and 3epsilon tests. Since each
original S_u is nonnegative, direct summation gives

    H_g(C) <= N*B_raw
              + eta*N(N-1)
              + epsilon*N(N-1)(N-2).                    (1)

No additional strict margin is needed. Empty and singleton codes are included
for the bounded-count certificate: all pair and triple factors vanish. A
fixed-N certificate still applies only at its specified N and radius cutoff.

## Comparison with the accepted repaired budget

The saved repair sets delta=epsilon*H+eta and uses

    B_M = B_raw+(M-1)delta,

where M=26 for the bounded certificate, M=N for a fixed-N certificate, and
H=T(1)=133770012193/4999950000. For the unchanged physical spherical energy
with positive common-g normalization W, replacing N*B_M by (1) increases its
proved floor by exactly

    G(N,M) = W*N*((M-1)delta-(N-1)eta
                            -epsilon*(N-1)(N-2))
           = W*[delta*N*(M-N)
                + epsilon*N(N-1)*(H-N+2)].              (2)

The count-repair constant is charged once in B_M and then its excess is
recovered by subtraction in (2). This is not another use of Q or an additional
spherical energy. Take the maximum between the old and new bounds from every
applicable saved certificate. Never sum gains from competing certificates.

For completeness, genuinely negative uniform repaired margins P'<=-p and
C'<=-c would additionally give the normalized gain

    N(N-1)p + N(N-1)(N-2)c/3.

That alternative uses repaired anchor positivity; it must not be added to
(2) without a separate accounting argument. The direct positive-error bound
(1) already replaces that repair accounting wholesale.

## Why this slack was not already consumed

The stage 2 driver (`../stage2_origin/original/verify.py`, lines 171 to 186;
the same code is in `../stage2_origin/run.py`) assigns the bounded proof surplus N*(Bbar-W*B_M), and the fixed-N
proof surplus N*(Bbar-W*B_N). It does not retain either term of (2). Its final
maximum is between those whole surpluses. Thus (2) is unspent in the accepted
profile. The original anchored PSD energies used in (1) are auxiliary proofs
of this same physical spherical energy, not extra physical components.

## Exact arithmetic artifact and scope

`run.py` in this folder runs the script `multiplicity/finite_count_gain.py`
on the stage 2 output (`python3 run.py --out DIR [--stage2 STAGE2_DIR]`).
The script reads the three accepted ANGULAR reports, checks the rational
repair and common-weight identities, and updates the saved 124 minima grouped
by (N, first-containing-radius index). All comparisons use Fraction.
Within one group the old and new applicable spherical bounds are constants,
so its minimum increases by the difference of those constants. This uses the
existing exact enumeration rather than rerunning 165503 count vectors.

This run is arithmetic from accepted reports, not a replay of their pair or
cyclic checks. The resulting JSON explicitly records that distinction.
