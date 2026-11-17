# Independent audit: 94D anisotropic parallelogram torus layers

**PASS**, as a periodic existence construction with occupied density at least **7.84434540405174e−21**. This is not an explicit-coordinate construction or a lattice-packing claim. No worldwide-priority guarantee is made.

## Existence and boundary issues

The Coxeter–Few–Rogers theorem applies to general Euclidean equal-ball coverings, hence to periodic lifted torus coverings. Its published simplex formulation was rechecked in G. Fejes Tóth, *Packing and covering in higher dimensions*, §2.2, https://arxiv.org/html/2202.11358 . Original theorem: https://doi.org/10.1112/S0025579300002059 . The regular-simplex cone integration in the inherited proof was reviewed; it uses n Dirichlet coordinates for the n-dimensional vertex cone, not n+1.

The exact first three moments of S=sum(T_i²) were independently expanded from Dirichlet moments. The cubic Taylor lower polynomial is globally valid because the fourth derivative of (1+S)^−46 is positive; the remainder has nonnegative sign on either side of the mean. The resulting rational covering lower constant equals the frozen one. All required forbidden-covering inequalities are strictly verified with an upper rational bound for pi.

If j closed torus balls covered the torus, their lifts would cover Euclidean space with density j V92 d^46/47^7, contradicting the strict lower covering bound. This reasoning does not require disjoint or injectively embedded torus balls. Their complement therefore contains a point at distance strictly greater than the selected radius. Choose a outside the one ball around0, then b outside the three balls around0,a,−a. This controls b,b−a,b+a simultaneously. The offsets0,a,b,a+b satisfy exactly the claimed two depth classes. Every difference between parity colors has coefficients0,±1 in a,b; all sign choices, including wraparound, are covered.

## Complete vertical and packing checks

The independent checker verifies the rational Gram and shows G−2I is positive definite. Enumeration of all nonzero pairs in[-3,3]² gives parity minima h1,h2,h2,4h1 for10,01,11,00 respectively. Any omitted integer pair has a coordinate of absolute value at least4, hence norm greater than32, exceeding all four recorded minima. Thus this finite enumeration is exhaustive, not a heuristic radius search. All16 ordered color cases meet squared distance16; same-color nonzero vertical differences have norm greater than16, while same-layer horizontal differences use the inherited parent minimum16.

The period L⊕2M has covolume4*47^7*sqrt(det M) and four translates. Their index cancels, giving covolume per point47^7*sqrt(det M). Radius2 in dimension94 yields center-density square2^188/(47^14 det M), checked exactly. Occupied-density rounding uses a rational LOWER Machin bound for pi and a squared rational comparison, not the producer's mpmath decimal.

## Parent and comparator provenance

The parent-source files parent_certificate.npz, unitary_sections.json and verify.py are byte-identical to the live public research/mw92_sections_20261115 package. The inherited verifier was freshly replayed:180 curve equations,16,110 pairwise polynomial-gcd heights, full integral generation, exact positive Gram and determinant47^14 all pass. Minimum16 remains the explicit published Schütt–Shioda Theorems6.23/13.35 and Section13.4.2 premise. We did not enumerate all parent vectors or claim to reprove that theorem. The archived book could not be fetched during this audit; source identification is retained from the existing public parent audit.

Live SpherePackings main was d2e55b3cd56c553beeb1d635aac5dbfd26b887dd. Its94D entry is approximately7.710930742674e−21. Comparing with its actual stored rational center-density square, rather than the rounded decimal, the new gain is strictly between1.73020% and1.73021%. The prior proof and exact row are preserved in baseline.json. The external table comparator MW92×A2 is explicitly a derived construction, not a surveyed global record; its exact comparison is also included in RESULT.json.

Files: check.py is an independently written rational checker; replay_parent.py reruns the inherited reviewed parent implementation. RESULT.json, PARENT_RESULT.json and PARENT_BINDING.json record the output and bindings. The frozen producer results.json hash is a7608134416525af4b8bf87c353c10af43dc754348b9541aa0b4034c1b2c62d4.
