# Search for further sphere-packing density records

## Outcome

**No new SOTA lower bound survived validation in this research run.** The
existing 36D construction remains the successful result. No searched 34D,
35D, 37–44D, 48D, or 56D candidate is promoted to a new record here.

The useful additions are reproducible searches, two exact exclusions of
specific construction families, and an exact certificate that the existing
36D lattice is a strict local density maximum among lattice packings. These
are mathematical and computational results, not new density records.

All comparisons use **center density**

```
delta = (minimum_squared_length/4)^(dimension/2) / covolume.
```

For an s-coset periodic packing, multiply the numerator by s. Packing fraction
is center density times the volume of the unit ball in that dimension. The
search objective is density, not kissing number.

## Baselines and provenance

The starting repository commit is
`db6b3f0e3c4c1ddd929a3dcb107d48a942af2498`.
Its 36D lattice has center density `16777216/3176523`, a factor
`139968/117649` over KP36, or an 18.970836981190% improvement.
The existing computer-assisted minimum certificate and its stated numerical
assumptions remain in `../../constructions/dim36/`.

For other dimensions the baseline is
[Henry Cohn's table](https://cohn.mit.edu/sphere-packing/), with dimensions
38–43 updated to the newer
[Dorofeev–Sun–Wang construction and certificates](https://github.com/wcgbg/sphere-packing).
The external checkout used here is pinned at
`a0bc7f51444e42c3cbae0d2a910039cd50013500`.
The cached Cohn table still contains older 38–43 values; those older values
are not the thresholds used in this search. This provenance check is not an
exhaustive world-priority guarantee.

## Results by route

| Route | Work completed | Outcome and scope |
|---|---|---|
| Existing 36D contact cone |94,050 exact antipodal contact witnesses; exact positive eutaxy identity; rank 666 modulo 65521 |Strict local lattice optimum, using the previously established minimum. No global optimality claim. |
|35D/34D sections of 36D | Short dual sections; metric cones; bounded searches for sections avoiding the norm 8 shell |No improvement. Six sampled35D cones closed numerically. A 23-contact parity obstruction is exact. |
|37D layers | Order3/order5 shifts, Voronoi vertex walks, joint37D Gram optimization |No improvement. Order-two extension is excluded exactly; the unrestricted searches are incomplete. |
|37D/35D/34D sections of 48D | Primitive normals of ranks11/13/14, including vectors of norms6/8/10/12 |Best tested normal determinants121500,391500,648675, all above record-beating thresholds. |
|38D/39D sections | Two pools with 300,457 and 600,454 exact vectors; primitive normal exchange |Best determinants54675 and 21870 reproduce the current records. |
|38–43D normal/cluster search |24,927 evaluated configurations; additional43D one-row screen with 10,703 distinct Gram matrices |No new candidate. Dual-cluster optimization is heuristic, with finite ray/point caps. |
|43D periodic deformation |946 metric +215 relative-shift variables;75,964 exact contacts available |A 12,000-contact numerical nonnegative-eutaxy solution closes the sampled first-order cone. |
|44D periodic deformation |990 metric +176 relative-shift variables; published P48n ambient;55,964 exact contacts available |A10,000-contact numerical nonnegative-eutaxy solution closes the sampled first-order cone. |
|48D modular seed | Six exact contact-null directions;546 vertices of a finite necessary-constraint relaxation |Best conditional density equals the existing 48D record. The discovery process does not itself certify the minimum. |
|56D quadratic pencil | One exact contact-null direction, integral Pell equivalence, two explicit short witnesses |Entire family has center density at most `(3/2)^28`, below the 56D benchmark. Exact, global within this family. |
|34–40D mixed codes | Exact outer-code enumerations and parity-column searches |The proposed34D rank-two fiber is impossible for the chosen code. Rank-three searches leave explicit bad supports. |

The 38–43D search count combines15,240 configurations in 41–43D and 9,687 in
38–39D. It is a count of evaluations, not a proof that these families have
been exhausted. The finite one-row screen tests all five replacements around
one43D normal basis against the saved600,454-vector pool, with the stated
heuristic cluster scoring.

## Exact certificates worth retaining

1. **Local optimality of the existing 36D lattice.**
   `contact_walk/verify36_local.py` verifies explicit minimum-vector norms,
   a positive rational identity for the inverse Gram matrix, and a full-rank
   modular minor. The proof and all positive weights are documented in
   `contact_walk/README.md`. This closes infinitesimal lattice deformation
   at the already successful endpoint, modulo scale and isometry.
2. **The 34D mixed-code parity obstruction.**
   The chosen `[17,9,7]` GF(4) code has a 2-(17,7,63) design of 408 weight-seven
   supports. Every admissible rank-two ternary parity coloring would force
   exactly 680/7 monochromatic coordinate pairs, a contradiction. This is
   checked independently in `mixed_codes/verify.py` and `sections/audit_n17.py`.
3. **The 56D entire-pencil bound.**
   `catalogue_screen/verify56_pencil.py` checks a compact integer certificate.
   Pell equivalence reduces every positive-definite parameter to 0≤t≤5, where
   two short vectors bound the density by 85222.69299239293. The reference is
   196213.21381197087. This exclusion does not depend on complete enumeration.
4. **The order-two37D extension barrier and finite section witnesses.**
   `gluing/verify_extension_thresholds.py`, `gluing/verify_normal_witnesses.py`,
   and `sections/verify_finite_results.py` check exact arithmetic, norm
   witnesses, determinants, saturation, and the small parity obstruction.

## Reproduction

The quick exact checks need Python, NumPy, and SymPy. Numerical discovery also
uses SciPy and fpylll/cysignals, plus a C++17 compiler. Use a single BLAS thread
per worker when running parallel searches.

From the repository root:

```bash
python research/density_breakthroughs_20261113/verify_exact.py
bash research/density_breakthroughs_20261113/build_tools.sh
```

The first command replays compact exact certificates and writes
`exact_verification_summary.json`; it does not repeat the exploratory searches
or label numerical eutaxy output exact. Each subdirectory README documents
its own search commands, caps, and assumptions.

Some discovery scripts use the external certificate repository. Recreate its
ignored checkout when needed:

```bash
R=research/density_breakthroughs_20261113
git clone https://github.com/wcgbg/sphere-packing "$R/seed_deformations/external_sphere_packing"
git -C "$R/seed_deformations/external_sphere_packing" checkout a0bc7f51444e42c3cbae0d2a910039cd50013500
```

Large contact pools, raw enumeration dumps, compiled executables, and the
external clone are excluded from git. Generator scripts, archived small source
inputs, compact witnesses, reports, and independent checks are retained.

## What would count as a next success

For a primitive rank-k normal lattice K in a published minimum-six unimodular
48D ambient lattice, a bare section has center density
`(3/2)^((48-k)/2)/sqrt(det K)`. Therefore an exactly checked determinant below
54675 at rank 10 improves38D, and below 21870 at rank 9 improves39D. Rank 11 needs
determinant below 102382.706… to improve37D. Any successful cluster must instead
pass the full exact expression

```
|S| * ((6-beta)/4)^((48-k)/2) / sqrt(det K),
```

with ambient membership, primitive saturation, distinct dual cosets, diameter
beta, and the current benchmark all checked. Those are candidate gates for
future work. No output in this checkpoint passes a new-record gate.
