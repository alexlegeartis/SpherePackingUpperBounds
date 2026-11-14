# Independent audit of the compact 36D and 56D claims

Both `../contact_walk/verify36_local.py` and `../catalogue_screen/verify56_pencil.py` were independently read and successfully rerun on the saved inputs. These are negative structural results for specified routes, not new packing-density records.

## 36D: exact extreme-form certificate

The reduced Gram satisfies Q=U Q_original U^T exactly, and detU=-1. Thus the existing minimum certificate applies under the advertised basis change. The94,050 displayed integer rows have exact Q-norm160. Their absolute coordinate values are at most14, making the intermediate weighted moment bound94,050·2581·14²=47,577,637,800 safely smaller than2^63; the verifier also checks the final multiplication bound.

The integer weights1326,2581,698,1326, assigned by the four exact deformation slopes, are strictly positive and satisfy

\[
\left(\sum_v c_vvv^T\right)Q=417590880I.
\]

Consequently Q^-1 is a positive linear combination of the contact outer products. The modular elimination proves rank666 in Sym36 over F65521, and therefore over the rationals. The prime check and modular multiplications are within exact signed64-bit arithmetic bounds.

These are eutaxy and perfection, respectively. With the prior minimum160 certificate, Voronoi's extreme-form criterion proves a strict local density maximum among lattices **modulo scaling and isometry**. Equivalently, in a fixed-minimum normalization the determinant is strictly locally minimized. Contact-shell completeness is unnecessary, because the explicit valid subset already supplies both identities.

This does not exclude a denser distant lattice, periodic or nonlattice packings, or constructions in other dimensions. Literal strictness without removing scale is not meaningful because packing density is scale invariant.

## 56D: global exclusion for the discovered one-parameter pencil

Write M(t)=G+tD. The exact verifier establishes G positive definite,

\[
\det G=29^{28},\quad H=G^{-1}D,\quad H^2=I/29,\quad \operatorname{tr}H=0.
\]

Since H is self-adjoint in the G inner product, its eigenvalues are real. They are +1/sqrt29 and -1/sqrt29, each with multiplicity28. Hence M(t) is positive definite exactly on |t|<sqrt29 and

\[
\det M(t)=(29-t^2)^{28}.
\]

The integral matrices U,R have determinant units and satisfy

\[
U^TM(t)U=\frac{27+5t}{2}\left(G+\frac{145+27t}{27+5t}D\right),
\qquad R^TM(t)R=M(-t).
\]

The scalar prefactor is positive throughout the positive-definite interval. Therefore these transformations preserve packing center density, allowing lattice basis change and positive rescaling.

For a direct finite-reduction proof, after reflection suppose5<t<sqrt29. Applying the inverse unit replaces t by

\[
g(t)=\frac{27t-145}{27-5t}.
\]

The denominator is positive. Moreover

\[
g(t)+5=\frac{2(t-5)}{27-5t}>0,
\qquad
g(t)-t=\frac{5(t^2-29)}{27-5t}<0.
\]

Thus either g(t) already lies in[-5,5], or the iteration remains positive and strictly decreases. It cannot remain above5 forever: its limit would be a fixed point of g in[5,t_initial], whereas the only fixed points are plus or minus sqrt29. It therefore enters[-5,5] after finitely many applications. Reflection yields0<=t<=5.

The two integer witnesses have M(t)-norms28 and32-4t. The actual lattice minimum is no greater than their minimum, so its center density is bounded by

\[
\left[\frac{\min(28,32-4t)^2}{16(29-t^2)}\right]^{14}.
\]

For0<=t<=1, the bracket is49/(29-t^2), increasing to7/4. For1<=t<=5 it is(8-t)^2/(29-t^2), whose derivative is

\[
\frac{2(8-t)(8t-29)}{(29-t^2)^2}.
\]

It decreases up to29/8 and increases afterward. Its maximum on that interval is attained at an endpoint and equals9/4, at t=5. Consequently the entire positive-definite pencil has center density at most

\[
(9/4)^{14}=(3/2)^{28}=85222.692992\ldots,
\]

which is strictly less than the supplied56D benchmark (5/2)^28/29^4=196213.213812...; the exact ratio is approximately0.4343371750.

The argument requires neither a complete minimum shell nor a claim that either witness always attains the true minimum. Upper witnesses are sufficient. Its scope is exactly this pencil and its integral equivalences; it is not a global upper bound on56D lattice packing.
