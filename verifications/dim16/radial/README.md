# Portable 16D radial certificate component

**Radial gates: PASS.** With spherical count U=7802, the exact outward-rounded
expression is strictly less than **0.024920988**. The unchanged U=7803 control
is strictly less than **0.024921003**. This directory proves the radial
component. The independently certified spherical count and full log-profile
cyclic geometry are separate theorem inputs, not silently assumed to have
been replayed here.

The U=7802 expression is

0.0249209879146585407778922410323795054326556519163892...

It is below the supplied comparator 0.02492121 by more than 2.22085e-7.
Both expressions and exact rational objectives are in
`expression_controls.json`. `RESULT.json` deliberately retains a component
scope instead of claiming this directory alone is the complete density proof.

## Run

Requirements: Python 3.12 and python-flint 0.9.0. No NumPy, SciPy, optimizer,
network, mutable-discovery import, geometry search, or external data is used.

```sh
python verify.py
```

The existing workspace environment can run it as:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  research/density_solid_cyclic_20261111/.venv/bin/python \
  research/density_log_count16_radial_certificate_20261111/verify.py
```

Do not run Python with `-O`. The verifier rejects disabled assertions. It
reconstructs the exact Laguerre relation, performs rigorous global root
isolation, proves the finite pair cover, checks negative controls, and writes
deterministic result files. The complete run takes about 0.2 seconds on the
current environment. The folder can be copied elsewhere and run unchanged.

## Exact function and one deterministic repair

The source degree-63 polynomial uses the normalized Laguerre basis
L_k^7(2X)/L_k^7(0). Its Fourier eigenvalue is (-1)^k. The source JSON numbers
are interpreted as exact binary64 rationals, and the verifier independently
reconstructs both monomial polynomials from the saved coefficient vector.
The log-profile coefficients are exact decimals, matching the separate
profile geometry and tail audit.

The source pair candidate is invalid: at exact
X=6858771052618467/562949953421312, its residual is strictly positive,
approximately 9.21252733225e-13. This failure is preserved and replayed as a
negative control, not hidden by sampled acceptance.

Only one repair is made:

1. Divide the source P, H, C, and lambda by the exact positive H(0).
2. Add epsilon=1/1000000 to Laguerre coefficient c1. Thus
   delta P=epsilon*(1-X/4), delta H=epsilon*(X/4-1).
3. Divide P, H, C, and lambda by 1-epsilon. The final H(0) is exactly 1.

The correction decreases P throughout the pair domain X>=2*pi and improves
H for X>=4. Its effect on the remaining Fourier segment is checked by exact
global root isolation. No LP or coefficient change is performed during
verification. The final exact polynomials and constants are in
`candidate_exact.json`.

## Global rays

Every isolated root has a rigorous Arb enclosure and the total multiplicity
equals 63. Real-root imaginary parts are required to be exactly zero.

- H has no roots on [0,infinity), and H(0)=1, so H is strictly positive there.
- P has no roots on [21/2,infinity), and P(21/2)<0, so P is strictly negative
  there. This is stronger than the requested physical tail starting at X=40.

All root boxes and exact anchor values are in `RAYS.json`.

## Complete pair inequality

Write v=X/(2*pi), h(v)=exp(-phi(v)). The proved inequality is

exp(-X) P(X) <= C/B*(2h-h^2) + lambda*1[1<=v<5000/4997].

The half-open threshold corresponds to squared distance
tau=10000/4997 and spherical cap 5003/10000. The count term is the ordinary
average-neighbor allowance, not another positive-semidefinite kernel. Its
objective cost is U*lambda without an edge-count factor of two.

On v in [1,7/4], every profile knot and the exact count threshold are split
before subdivision. On each analytic piece, phi(v)=a*v+b*log(v)+c.
The verifier multiplies the residual by exp(X), evaluates its midpoint and
first derivative, encloses the complete second derivative with Arb, and
uses Taylor's theorem to bound the entire interval. Polynomial derivative
evaluation uses centered polynomial composition to control cancellation.

All 36 accepted rational intervals have strictly negative upper bounds.
They form an exact contiguous partition of width 3/4; maximum depth is 6.
The counted left limit at 5000/4997 and uncounted right limit, including the
threshold itself, are both proved. No floating-point sign acceptance occurs.

For v>=7/4, X>=7*pi/2>21/2 and P is already negative. Nonnegative profile
coefficients imply 0<h<=1, hence all allowances are nonnegative. This covers
the entire remaining pair ray, including all later kinks and the requested
finite endpoint X=40. The saved partition is in `partition.json`.

## Negative controls and hashes

`negative_controls.json` records:

- the exact positive residual of the unrepaired source;
- rejection of a cover with the first interval deleted;
- rejection of incorrect closed-boundary count semantics.

`SHA256SUMS` and `manifest.json` identify the frozen files. The saved numerical
source is included only to bind the exact rational polynomials to their
Laguerre construction; numerical success flags are not used for acceptance.

## External theorem components

The integration must bind these separately certified components:

- Spherical counting gate: `density_kissing16_independent_cyclic_20261111`,
  exact cap 5003/10000, exact real objective about 7802.2151671287232,
  integer count at most 7802.
- Fixed angular and log-profile geometry:
  `density_log_hinge_geometry16_20261111`, including compact geometry,
  all-scale reduction, and infinite cyclic tail.

Their report hashes and this component's matching profile/angular hashes
are recorded in `external_dependencies.json`. Those report snapshots are
provenance, not a replacement for independently replaying their proofs.
