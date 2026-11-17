# Independent audit of repaired degree29 general50 bound

PASS: scalar, PSD, corner signs, normalization and complete512bit trace replay all pass. All32905nodes/16453leaves were independently replayed with no uncovered boxes. Safe unrestricted density upper bound <1.99548924e-7. See main_independent_result.json.

## Exact input bindings
Geometry frozen.json and scalar kernel.json are byte-identical, SHA2562d4ffd3c113149707bf2b55aff7986c666c75666c62bd5ad9b797ac005a849c2. check_bindings.py verifies this, radial polynomial identity, every annular endpoint/adjacency and optimized cap use for each positive height. The safe bound is strict: exact final rational density <199548924/10^15.

## Kernel/theorem audit
For degree k the normalized Gegenbauer addition kernel is PSD in50dimensions. Each positive c_k and exactPSD2x2B_k preserves this property; the origin block [[s,z^t],[z,B0]] PSD closes the specialorigin. c0=1 and q=sumc normalize the origin/offdiagonal consistently. Kdiag=(H^tWH)/q, mixed2K(0,v)=2z^tH/q. There is no lattice, mod2, antipodal or doubled diagonal assumption.

The compact sorted squared-side parametrization maps the entire cube to1<=a<=b<=c<=19/13. The expression checked by the gate is the cyclic sum divided by h(a)h(b) in the interior. Multiplying back gives a continuous identity through zero-factor boundaries, so the boundary is covered without dividing by zero. All three sides form a triangle since sqrt(19/13)<2. If exactly the longest side crosses support, the surviving angle lies in[-1,1/2]. The three corner polynomials bound the entire bilinear rectangle in h(a),h(b)∈[0,1]. If two sides cross, the cyclic sum vanishes. Corner-polynomial global extrema are enclosed by exact rational polynomials plus512bitArb enclosures of all complex derivative roots and both endpoints.

Direct periodic PSD summation plus cyclic triangle grouping yields numberdensity <= f(0)+s+E. CompactK enters finite local sums and needs no Fourier positivity; f remains a Schwartz Gaussian times polynomial and its Fourier transform is certifiednonnegative. This establishes an unrestricted packing bound, using the standard reduction toperiodic packing densities. No new restriction is being imported from lattice-only certificates.

## Interval implementation review
jet.py's product,reciprocal,square-root,and integer-power derivative/Hessian formulas were reviewed. pe2 is a centered polynomial Horner enclosure with simultaneously differentiated first/second derivatives, composed via chainrule. The natural interval, first-order mean-value and second-order Taylor bounds used by the geometry gate all enclose the value over each dyadicbox. Floating scores select partitions only; they do not certify signs. The replay independently reads each tracebranch, splits its rationalbox, verifies all leaves at512bits (producer384), and requires a fully exhausted partition with nodes=2leaves−1.

## Scalar proof replay and review
Fresh exact factorization and closed-form Fourier involution pass for degree71; Q(0)=1. Negative-linear factors1+x/a havea>0; positive-double factors are squares; conjugate-pair factors are globallypositive. Hence Q(x)>=0 forx>=0.1024bitArb accountsforall71spatial roots, allpossible realroots below180; exactP(180)<0 implies fullnegative tail.

Fresh768bit radialreplay proves633annuli covering[22,180],883partitionnodes. Its polynomial-shift expansion,20term exponential Taylor remainder and concavequadratic bound were reviewed. Nonnegative maxheight perannulus upper-bounds residual. Kernel support boundary is an explicitedge. All positiveheight caps have exact Delsartepolynomial reconstructions, with sign supplied by(t−s)(t+1)^bK(t)^2 and nonnegative Gegenbauercoefficients. Zeroheightrows contribute0, so analytic fallbackcaps are not needed. For squared radii ina narrowannulus [l,r] withr−l<=d, separated points normalizedtosphere haveinnerproducts <=1−d/(2r); monotonicity ineach radialvariable followsfrom|a−b|<=d. Thus caps apply to unrestricted neighboringcenters.

Optimized errorbudget≈0.00003362046733717; exactdensity≈1.9954892364793078e-7. Contactx=d=22 correspondsdistance²22/pi; volume-normalization cancels pi and equals(22/4)^25/25!. This is the scalar verifier's exactfactor.

## Audit scope
The numerical discovery candidate is not used. Scalar verifiers and intervalJet/gate were reviewed and freshly executed, not wholly reimplemented; inputbinding/normalization checks are independently written. No publication performed by this audit.
