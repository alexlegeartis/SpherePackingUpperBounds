# Anisotropic parallelogram torus layers in94D

Candidate awaiting independent audit. This is an existence construction, not explicit coordinates for horizontal offsets.

Use the published min16, covolume47^7 MW92 parent L and the same inherited Coxeter–Few–Rogers covering lower bound C as the existing cubic-covering package. Let d_j be rational depths satisfying j V92 d_j^46 < C47^7. Such j equal torus balls cannot cover the torus: their lifts would contradict the covering theorem.

Choose a in R92/L with distance²(a,L)>d_1. Next choose b outside the three radius-sqrt(d_3) torus balls centered0,a,−a. Then b,b−a,b+a each have distance² to L>d_3. Consequently the four offsets0,a,b,a+b have horizontal squared separation at least d_1 for the matching pairs differing by a, and at least d_3 for all other pairs.

Set h1=16−d1 and h2=16−d3. Let the vertical rank2 lattice M have Gram

    [[h1,−h1/2], [−h1/2,h2]].

Color its points by coordinates modulo2 and assign the offsets0,a,b,a+b to00,10,01,11. The period is L orthogonally summed with2M, and there are4 translates per period. For integer coordinates x,y, the vertical norm is

    h1(x−y/2)^2+(h2−h1/4)y^2.

When y is odd this is at least h2, attained for y=±1 and nearest integral x. If y is even and x odd, the norm is at least h1: for y=0 immediate; for |y|≥2 it is at least4h2−h1>h1. Nonzero vectors with x,y even have norm at least4h1>16, since the lattice minimum is h1. Thus every horizontal/vertical difference has squared norm at least16. Same-layer parent vectors attain16.

Covolume per point is47^7 sqrt(detM); the4 translates cancel the index4. Center density squared in94D is

    (2^184/47^14) * 16/detM.

The exact rational checker test.py derives C from exact Dirichlet moments, brackets pi above using alternating Machin series, verifies the strict forbidden-covering inequalities, and supplies the exact rational Gram and density. The decimal occupied density is approximately7.8443454040517406e−21. It is not yet an outward-rounded registered bound. Gain over the existing uniform4-color Gram is approximately1.730201766%.

The novelty relative to the earlier uniform coloring calculation is the structured offset parallelogram and nonuniform parity metric; the covering theorem and parent are inherited inputs. No worldwide-priority claim is made here.
