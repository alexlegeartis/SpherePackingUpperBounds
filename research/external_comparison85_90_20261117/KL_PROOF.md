# Stronger external-derived MW92 section comparisons

These are new numerical evaluations of **pre-existing methods**, not new sphere-packing methods or dimension-specific records claimed to be printed in the cited papers. They supersede the weaker Rogers-derived comparisons in89D and90D. No global optimization claim is made.

## Published ingredients

[Cohn–Zhao, Sphere packing bounds via spherical codes (2014), Proposition2.1](https://arxiv.org/html/1212.5966) states that for pi/3≤theta≤pi, the packing density in R^n is at most sin(theta/2)^n A(n,theta). Their Theorem3.1 states the spherical polynomial LP bound. [Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini (2020), §5](https://doi.org/10.1007/JHEP12(2020)066) discusses the Levenshtein spherical bounds and this transfer. The MW92 parent minimum16 and determinant47^14 come from [Schütt–Shioda, Theorem13.35](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf), already used and attributed in the packing package.

## Exact finite certificates

For each n=90,91,92, check_kl.py constructs a degree11 polynomial f(x)=(x−s)K(x)^2 using a Jacobi reproducing kernel, at the rational s recorded in kl_report.json. It verifies the normalized Gegenbauer expansion has nonnegative coefficients and positive constant coefficient a0. Since s lies in(−1,1/2), the transfer theorem applies. The square factor proves f≤0 on[−1,s] without numerical root assumptions. Thus A(n,theta)≤f(1)/a0 and the squared density upper bound is

    U_n² = (f(1)/a0)² ((1−s)/2)^n.

The independent verifier does not reuse the Jacobi generation or norm formula. It reconstructs normalized Gegenbauer polynomials by their three-term recurrence from the saved coefficients, divides by(x−s), and factors the quotient exactly to prove that it is a positive scalar times a square. It separately checks all coefficients, rational bounds and integer powers.

The Hermite constant satisfies gamma_n≤4(U_n/V_n)^(2/n). Rational Machin bounds on pi give certified values

| n | Rational upper bound on gamma_n |
|---:|---:|
|90|11996914884433/1000000000000|
|91|3026811744289/250000000000|
|92|488702762957/40000000000|

The original verifier uses100/101 alternating arctangent terms; the independent verifier uses120/121. Both verify the final power inequalities exactly.

## From Hermite bounds to packing lower bounds

For an n-dimensional lattice of squared minimum m and covolume d, a primitive shortest dual vector v has ||v||²≤gamma_n d^(−2/n). The perpendicular section has covolume d||v|| and minimum at least m. If D is the original center-density square, the section center-density square is at least D^((n−1)/n)/gamma_n. The dual vector is primitive because a nontrivial integer multiple cannot be shortest.

Starting with D=2^184/47^14 for MW92, apply this formula in92,91,90 dimensions. Each step rounds DOWN to an integer q, with the exact check

    q^n gamma_n^n ≤ D^(n−1) < (q+1)^n gamma_n^n.

The resulting certified center-density-square floors are:

| Section dimension | Floor |
|---:|---:|
|91|3511996470050690436868669299996|
|90|133919321295115171411465916890|
|89|5298367406746324075995062827|

Multiplying their square roots by the unit-ball volumes gives approximately7.200357004613e−20 in90D and5.435520207304e−20 in89D. These are achieved lower guarantees (existence of sections), not Rogers-style upper endpoints. Our comparison gains use these exact floors, not the rounded decimals. The fixed rational angles and degree are not asserted globally optimal. Stronger published numerical Euclidean LP bounds may further strengthen this external-derived comparison.

## Replay

    python research/external_comparison85_90_20261117/check_kl.py
    python research/external_comparison85_90_20261117/independent_kl.py
    python research/external_comparison85_90_20261117/verify.py

Requires Python and SymPy. All three checks passed on17 November2026. No packing certificate or claimed density of our constructions was changed.

## Stronger numerical data are not silently treated as certificates

The authors' [MIT dataset](https://hdl.handle.net/1721.1/125646), dataminusone.txt, includes numerically optimized scaling gaps for90–92D. The highest-N values are approximately5.737208410825264,5.792122865905022 and5.847008716053789. Under the paper's normalization these suggest Hermite bounds twice those numbers, stronger than the fixed-angle certificates above. They are numerical spectra, not the verified auxiliary-function positivity certificates needed for an exact section guarantee. The paper's §2.4 distinguishes that positivity step. We do not claim to have certified these stronger numerical values or optimized the external baseline. Their existence is disclosed rather than calling the weaker analytic evaluation globally strongest.
