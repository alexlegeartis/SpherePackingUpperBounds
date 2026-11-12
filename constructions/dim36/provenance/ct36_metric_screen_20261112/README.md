# 36D sparse-shell metric screen

No new packing record. Source: Nebe–Sloane catalogue CTturyn2.html, Delta2, dimension36, determinant3^18, minimum6, kissing576.

## Independent checks and scope

Parsed source Gram is symmetric integral even and its exact determinant is387420489=3^18. A unimodular Gram-LLL change of basis produced H. Completed floating Cholesky enumeration through norm6 found576 vectors; every returned norm was checked exactly to equal6. This matches the catalogue, but floating branch completion is NOT an independent exact proof of absence of shorter vectors. The family exclusion below requires only the explicit vectors, not that absence claim.

## Three-block reduction

Let V be the576 checked norm6 rows and A=V^T V H. Exact matrix arithmetic checks A^3=288A^2-18432A; A has eigenvalues0,96,192 each with multiplicity12. Its H-orthogonal projectors P0,P1,P2 define the three-positive-parameter Gram
H(a,b,c)=H(aP0+bP1+cP2).
Thus det(H(a,b,c))/det(H)=(abc)^12.

At minimum8, center density relative to the current36D record is3/(abc)^6. Beating the record therefore requires(abc)^12<9. Individual norm6 vectors have block squared lengths(0,2,4).

## A tempting candidate and its exact failure

The original norm6/8 shell constraints admitted(a,b,c)=(1/3,1,3), an even integral Gram with unchanged determinant. A separate global floating enumeration below8 found576 vectors of exact candidate norm6. Their original norm is14 and their block lengths are(12,2,0). They give the essential missing inequality12a+2b>=8. These explicit witnesses rigorously reject the candidate independently of enumeration completeness.

## Exact exclusion of the whole three-scalar family

Twenty-one explicit integer vectors supply rational constraints qa*a+qb*b+qc*c>=8. Exact checks recompute their block lengths from the integral matrices. In particular, vectors of block lengths(24,0,0),(0,8,0),(0,0,8) force a>=1/3,b>=1,c>=1. The feasible polyhedron has four vertices, enumerated with rational Cramer determinants. At every vertex abc>=5/4. The recession cone is the nonnegative orthant, and log(abc) is concave and increasing in each coordinate, so its minimum is attained at a vertex. Hence abc>=5/4 throughout the feasible set.

Consequently every metric in this three-block family has density at most3/(5/4)^6=12288/15625=78.6432% of the current36D record. This is a necessary-constraints upper bound, so additional unknown short vectors cannot invalidate the exclusion. Homogeneous rescaling covers all minimum normalizations, not just minimum8.

This does not exclude general666-parameter Gram deformations, other polarizations, or nonlinear packings. No general metric optimizer was run after this cheap structured family failed.

## Reproduction

exact_screen.py verifies explicit-vector constraints and enumerates rational vertices; it currently reads saved data under/tmp/ct36. Replace that prefix with this directory to run from the archive. All required witness arrays and source Gram are included. The bulky exploratory norm8 shell is unnecessary for the exact exclusion.
