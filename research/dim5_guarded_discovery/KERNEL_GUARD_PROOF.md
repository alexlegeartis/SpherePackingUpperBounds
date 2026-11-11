# A strict PSD kernel correction with zero pair cost

This is a verification aid, not a new density hierarchy. It can repair tiny
far/cyclic sign violations while preserving the radial function and complete
pair inequality. The isolated-origin cost may make large repairs useless.

For support radius R=3/2 and center separation one, every surviving kernel
angle is at most 7/9. Indeed, (r^2+s^2-1)/(2rs) is separately convex for
1<=r,s<=R, so its maximum is attained at a corner, where the candidates are
1/2, R/2=3/4 and 1-1/(2R^2)=7/9.

Let P(t)=sum_l p_l G_l(t), with p_0=1, p_l>=0, and P(t)<0 throughout
[-1,7/9]. Such a polynomial is obtained by a finite LP and then checked with
exact rational Sturm arithmetic. Set F=P-1/2, so all its Gegenbauer
coefficients c_l remain nonnegative, c_0=1/2, and F(t)<-1/2 on that interval.
Write D=F(1)=sum c_l.

Use annulus-constant radial features. Set the nonzero block in harmonic zero
to the 2-by-2 matrix

    [[D^2/(4c_0), -D/2], [-D/2, c_0]],

where the first coordinate is the isolated origin. This matrix is PSD with
zero determinant. In each higher harmonic degree, add c_l only in the constant
radial feature. The addition theorem therefore proves the complete guard
kernel positive semidefinite.

Its pair correction is exactly

    2 K_guard(x,0)+K_guard(x,x) = -D+D = 0

on the closed annulus, and zero outside it. Every surviving off-diagonal
kernel value is less than -1/2. Thus its compact cyclic sum is less than
-3/2, and its one-long-edge value is less than -1/2. Configurations with two
long edges have no surviving terms.

For any nonnegative rational epsilon, adding epsilon times this kernel
preserves PSD and the pair inequality and decreases every nonzero geometric
expression. The isolated-origin cost before possible further Schur
saturation is epsilon D^2/(4c_0). Exact full verification must still be run on
the final source; the numerical choice of epsilon is not proof evidence.

The current implementation only discovers exact guard polynomials for
support radius3/2 and specified dimension. It does not automatically repair
any witness and does not prove a density bound by itself. No published count
input is needed for this correction lemma.
