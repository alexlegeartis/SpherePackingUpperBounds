#78D three-fiber projection with a162/41 complex-multiplication normal pair

Status: primary exact replay passed; independent audit pending. This is a periodic packing, not a claimed78D lattice or an upper bound.

The parent is the rank80 MW lattice of y²=x³+1+t42 overF41², with minimum14 and Gram determinant41^12. Its190 supplied polynomial sections, exact height Gram, integral generation, determinant and positive definiteness are replayed by parent/verify_parent.py. The height lower bound follows fromχ=7,42 irreducible typeII finite fibers and smooth infinity, as in Schütt–Shioda's Mordell–Weil Lattices (2017 author draft), Chapter13: https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf . This published height theorem is an external mathematical input.

The two exact integer functionals B have normal Gram

    K = [[162,-81],[-81,162]]/41.

Smith form1,1 proves the parent-to-label map is ontoZ². Use labels S={(0,0),(1,0),(0,1)} and project each labeled fiber onto the common orthogonal kernel. Every nonzero pairwise label difference has squared normal projection82/243. Consequently every projected distance is at least

    14−82/243 =3320/243.

This bound is attained: source parent minimum vector2 has label(1,0). Its projected norm is exactly3320/243 and joins the zero fiber to the first nonzero fiber. No shortest-vector enumeration is used or assumed complete.

The saturated kernel period has Gram determinant

    det(G)det(K)=3^9*41^10=264198203201729708883.

freeze78.py constructs a unimodular HNF transform, the full78×78 period Gram, three rational shifts in its basis, and the attaining vector. Distinct fibers cannot coincide after projection because their projected mutual distances are positive. Thus there are exactly three translates per period.

The center-density square is exactly

    9*((3320/243)/4)^78/(3^9*41^10)
    = (830/243)^78/(3^7*41^10).

A rational Machinπ enclosure proves ordinary density at least1.41519803035180e−16.

## Comparisons

The live repository README, retrieved17November2026 after its09:15:26UTC modification, still displayed78D1.366254130898e−16, the earlier MW80 antipode. Its exact center-density square is pinned in candidate78.json from research/range48_96_20261115/results.json. The new-to-prior ratio is1.0358234228511676, about3.58234% improvement.

The source-attributed external derived comparison remains L80 perpendicular to2A2, from Bachoc–Nebe (1998): https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/extremm22.pdf . Its exact center-density square is2^74/3, ordinary density≈9.523775961404e−17. The new ratio is1.4859631684817027, about48.59632%. This is a derived comparison from the cited lattice and its previously audited normal section, not a78D record printed in that paper. No exhaustive worldwide-priority claim is made.

## Replay

Run PYTHONPATH=/tmp/packing-deps python freeze78.py in the present environment, or install numpy and python-flint and run python freeze78.py. It replays the parent and all geometry, determinant, minimum-attainment and rigorous-density checks. build78.py records discovery dependencies and is not needed for the frozen proof replay.
