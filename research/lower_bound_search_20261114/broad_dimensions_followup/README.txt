Broad-dimension lower-bound follow-up, 14 November 2026

No new sphere-packing density record is asserted. The certified 36D construction
remains the repository's lower-bound improvement. This package explores other
dimensions and supplies exact exclusions and a verified local deformation
obstruction. The earlier general self-orthogonal minimum-12 problem remains open.
All densities in this package are CENTER densities; comparisons use their exact
rational squares.

The main outcomes are:

* In dimensions 19, 20, 21, 22, 44, 45 and 46, an exhaustive finite search over
  minimum-basis section Grams and arbitrary finite dual clusters matches the
  existing antipode record. Every squared cluster diameter 0 <= beta < mu is
  covered. This includes clusters larger than simplices.
* In 43D the analogous rank-five search matches the existing record for
  beta <= 13/20. It checks 51,647 relaxed Gram types and solves 1,571,848 clique
  problems, with no node-capped cases left unresolved.
* In 22D and 46D, a separate search also allows longer rank-two bases: all even
  reduced Grams [a,b;b,c] with mu <= a <= c <= 64 and 0 <= b <= a/2. All cluster
  diameters and sizes are covered. The 6,448 and 6,355 forms give no improvement.
* For the five-translate 44D record, 13,000 exact contacts give a full-rank,
  strictly positive stress certificate. Metric changes and independent motions
  of all five translates have no first-order density improvement.
* One million order-three and one million order-four translates of the 36D
  record were screened for a 37D extension. All were rejected. The 92,275
  promoted samples have archived exact nearest-point upper witnesses.
* A 37D full-metric LP found a direction that preserves the saved contacts to
  first order. New short vectors invalidate the tested finite steps along it.
* Thirty-five catalog lattices in 28D, 30D, 31D and 32D were screened. The Q30
  and Q31 minimum-shell eigenspace families give zero first-order descent.
  These spectral screens and their floating shell completions are exploratory.

Exact anticode search

Take a primitive section K of rank k in an integral unimodular N-dimensional
lattice of minimum mu. Assume K has a lattice basis consisting of vectors of
squared norm mu. This is a condition on the saturated section's basis; a
sublattice spanned by independent shortest vectors need not be saturated.
For S in K* with s points and squared diameter beta < mu, the antipode score is

    delta^2 = s^2 / det(K) * ((mu - beta)/4)^(N-k).

The searched cases are (N,mu,k) = (24,4,2..5) and (48,6,2..5).
Every eligible Gram has diagonal mu and integer off-diagonal entries between
-mu/2 and mu/2; otherwise a difference of two basis vectors would be shorter
than mu. Enumerate these matrices modulo signed coordinate permutations. Reject
matrices with a nonpositive leading principal minor or an explicitly short
vector in {-1,0,1}^k. The survivors are a SUPERSET of eligible sections, so
discarding impossible matrices is unnecessary for an upper exclusion unless
one of them produces an apparent gain.

The rank-five relaxation does admit three forbidden Gram types in each ambient
minimum. They contain a vector with coefficient 2 and squared norm 2, missed by
the first finite filter. Any apparent gain is now subjected to a complete exact
minimum check before it can update the score. For an even integral Gram, a
short vector has norm at most mu-2, and Cauchy-Schwarz gives

    |v_i| <= floor(sqrt((mu-2) * (K^-1)_ii)).

The resulting finite integer box is checked exactly. The rejected witnesses
are archived in rank5_results.json and rank5_mu6_results.json. The original
relaxed gains are not constructions or density records.

Translate a cluster so that 0 is one of its points. A dual point has integer
pairing coordinates h, squared norm h K^-1 h, and

    |h_i| <= sqrt(beta * mu) < mu.

For the full search it therefore suffices to enumerate h_i = 1-mu,...,mu-1.
For each possible diameter, make the graph joining two eligible points when
their squared difference norm is at most beta. A maximum clique, together with
0, is an optimal cluster. All graph arithmetic and density comparisons are
integer or rational. Before solving a graph, use |S-S| >= 2|S|-1: if a dual
ball contains j nonzero points, s <= j/2+1. This prunes large-diameter cases.
The graph search has an explicit node cap and records capped cases as UNKNOWN;
the archived exact exclusions have no UNKNOWN cases.

The all-diameter searches solve 3, 35, 712, 69,400, 5, 87 and 10,385 graphs in
22D, 21D, 20D, 19D, 46D, 45D and 44D respectively. The rank-five Gram generator
has 6,595 signed-permutation orbits for mu=4 and 160,419 for mu=6 before the
finite validity filter. An independent Burnside count checks these totals.

These exclusions do not cover sections without a basis of minimum-length
vectors, different ambient minima, or general packings. The 43D diameter range
is restricted. No global optimality claim for sphere packing follows.

Rank-one antipodes, including longer sections

For completeness, rank one has a short analytic optimum. A section generated
by a primitive vector of squared norm a has dual spacing 1/sqrt(a). An s-point
cluster, s>=2, has beta >= (s-1)^2/a. Put z=(s-1)^2/a. Its density is at most

    s/(s-1) * sqrt(z) * ((mu-z)/4)^((N-1)/2).

The maximum over 0<z<mu occurs at z=mu/N. Since s/(s-1) decreases for s>=2,
s=2 is best. The required a=N/mu is 6 in the Leech case and 8 in the P48p case;
these sections exist and give the known 23D and 47D records. The s=1 optimum,
at a=mu, is slightly smaller in both cases (verified by exact squared-density
comparison). Thus more points or longer rank-one sections cannot improve these
two records within the same ambient construction.

44D positive-stress certificate

antipode44_0.npz specifies a primitive rank-four section in an archived P48p
basis, its 44D complement Gram Q, and five shifts with denominator 55.
The ambient minimum six is a published input, not proved anew here. The
complement has determinant 605; the cluster diameter is 16/55 and the packing
minimum is m=314/55. Each stored contact is v=w/55, for integer w, with

    w Q w = 17270 = 55^2 * m.

The verifier reconstructs the primitive complement and shifts with integer
arithmetic, checks unimodularity, positive definiteness and every contact's
coset congruence and exact norm. No random sampler or LP solution is trusted.

For nonnegative contact stresses lambda, the necessary equations are

    sum_i lambda_i w_i w_i^T = 3025 Q^-1,
    sum_i lambda_i sign(site,i) w_i = 0  for sites 1,...,4.

Fix site 0 to remove global translation. The 990 symmetric-metric equations
and 176 translation equations form C lambda = t, with 1,166 rows. Both C and t
are integer here. A strictly positive exact solution, and full row rank, give
the obstruction: every nonzero first-order deformation preserving the active
distances increases the covolume to first order.

The certificate establishes existence of an exact rational positive stress
without printing enormous rational weights. Let D be a diagonal power-of-two
row scaling, A=DC, b=Dt, G=A A^T, and let H and w0 be the stored/computed binary
rational approximations. Define E=I-GH, R=A^T H and r=b-Aw0. Then

    lambda = w0 + R (I-E)^-1 r

satisfies A lambda=b exactly. If ||E||inf<1, the correction has norm at most

    ||R||inf * ||r||inf / (1-||E||inf).

The verifier bounds every matrix product using the standard binary64 dot
product error bound, with epsilon=2^-52, twice the unit roundoff, and additional
outward inflation. The integer Gram product is exact: its largest row square
sum is 664752643826350 < 2^53. Cauchy-Schwarz bounds every dot product's sum of
absolute integer terms by this number, so every product and partial sum is
exactly representable in binary64. Power-of-two scaling is exact as well.
An approximate inverse is accepted only after its residual is verified.

The archived/replayed bounds are conservatively summarized by

    ||E||inf < 0.001422,
    ||R||inf < 649,
    ||r||inf < 2.577e-12,
    correction < 1.675e-9,
    min(w0) > 3.98e-5.

Consequently the exact stress defined above is strictly positive.
The same residual bound proves G invertible and hence all 1,166 equations
independent. This is a verified obstruction for this five-site periodic model.
It says nothing about changing the number of translates or moving to a remote
packing.

37D searches

Stacking the 36D record L with generator (x,h) has center density
delta36*sqrt(2)/h at squared minimum eight. Beating 4*sqrt(2) requires

    dist(x,L)^2 > 8 - delta36^2/16 = 6.256524693309225...

The order-three and order-four screens used caps for discovery only. Promoted
samples were rejected by explicit lattice points: the numerator norms are at
most 1120 and 2000 in denominators 3 and 4 respectively. Dividing by 20*q^2
gives squared-distance upper bounds 56/9 and 25/4, both below the threshold.
The exact auditor checks all 47,444 and 44,831 archived promoted witnesses.
This is not a covering-radius upper bound or an exhaustive search of the torus.

The full 37D metric direction is also only exploratory. One tested step has
conditional center density 5.8630 at the assumed minimum eight, but an explicit
layer-one vector has squared norm about 6.09868. The conditional gain is invalid.

Reproduction

Required: Python, numpy, scipy, python-flint, and a C++17 compiler. fpylll and
cysignals are additionally needed for the optional catalog screen.

    OPENBLAS_NUM_THREADS=1 python verify_all.py

This uses temporary files, independently tests clique correctness on small
graphs, checks the Burnside counts, regenerates both rank-five Gram archives,
repeats the short exact searches, checks the 44D stress certificate, and audits
the 92,275 promoted deep-hole witnesses.

    OPENBLAS_NUM_THREADS=1 python verify_all.py --full

This also repeats the longer 19D all-diameter and 43D bounded-diameter searches.
The recorded 43D run took about eleven CPU minutes; runtime depends on hardware.

Optional exploratory replays:

    OPENBLAS_NUM_THREADS=1 python catalog_screen.py
    OPENBLAS_NUM_THREADS=1 python spectral_screen.py
    OPENBLAS_NUM_THREADS=1 python metric37.py
    OPENBLAS_NUM_THREADS=1 python ray37.py

The catalog and metric scripts generate files in work/. Floating enumeration,
spectral LP outcomes and CVP output are not minimum certificates.

Sources

Chen, Hu, Li, Wang, Wu, New Sphere Packings from the Antipode Construction:
https://arxiv.org/abs/2505.02394

Dorofeev, Sun, Wang, Optimal Extensions of Cross-Sections: Sphere Packings in
Dimensions 38 to 43, version 4:
https://arxiv.org/abs/2607.20359v4

Ambient P48p certificate source, pinned to the same revision as the earlier
research package:
https://github.com/wcgbg/sphere-packing/tree/a0bc7f51444e42c3cbae0d2a910039cd50013500

Nebe--Sloane catalog, archived HTML seed pages in catalog_sources/:
https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/

The original 36D lattice and its complete minimum proof are in
../../../constructions/dim36/.
