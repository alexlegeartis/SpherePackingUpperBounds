# 94–95D: anisotropic torus-offset flags with parity-weighted root chains

Verified occupied-density lower bounds:

| Dimension | Safe lower bound | Gain over preceding public certificate | Method |
|---:|---:|---:|---|
|94|≥7.84434540405174e−21|1.73020–1.73021%|Parallelogram torus offsets with anisotropic planar parity layers|
|95|≥2.06770138493108e−21|1.2018716893–1.2018716894%|Recursive binary torus-offset flags with weighted A3 parity layers|

These are **periodic existence constructions**, not explicit-coordinate or lattice-packing claims. Coxeter–Few–Rogers covering bounds guarantee horizontal offsets; vertical Gram matrices and all separation inequalities are exact. No new96D result or exhaustive worldwide-priority assertion is made.

## Reproduction

From any directory run `python verify.py`. Requirements: Python 3, NumPy and python-flint. The new94/95 arithmetic checkers use only the standard library; the inherited parent replay additionally uses NumPy/python-flint. The verifier checks a full manifest, both rational construction proofs and outward rounding, and freshly replays the public parent-source checker. It uses the existing ../mw92_sections_20261115/{verify.py,parent_certificate.npz,unitary_sections.json}; their hashes are bound by the manifest. No network or solver is needed.

- [94D independent audit](REPORT.md), [95D independent audit](REPORT95.md)
- [Parallelogram proof](input/PROOF.md), [recursive flag proof](input/FLAG95_PROOF.md)
- [94D exact result](RESULT.json), [95D exact result](RESULT95.json)
- [Parent-source replay](PARENT_RESULT.json)

The frozen producer scripts in input/ are retained for provenance, not used by verify.py; their decimal displays use mpmath, while the independent checks are exact fractions. The original proposal wording “pending audit” is historical and superseded by the independent reports here.

## Sources and scope

The MW92 parent minimum16 and determinant47^14 are inherited from Schütt–Shioda, Theorem 13.35 and Section 13.4.2. Its explicit section/Gram identification is [replayed separately](../mw92_sections_20261115/verify.py). The determinant, full integral generation,180 curve equations and16,110 pairwise polynomial-gcd heights are freshly checked; the published minimum theorem is not reproved by enumeration.

The covering theorem is [Coxeter–Few–Rogers 1959](https://doi.org/10.1112/S0025579300002059), stated in the [author survey §2.2](https://arxiv.org/html/2202.11358). The exact cubic Taylor/Dirichlet-moment covering constant is inherited from [the preceding proof](../adaptive_layers_20261115/COVERING_PROOF.md) and independently re-derived by check.py. Strict forbidden-covering inequalities handle overlaps and boundary points without requiring torus balls to embed injectively.

The prior public94/95 certificates remain untouched; baseline.json and baseline95.json pin their exact values and commit. The main repository table retains its external **derived** MW92×A2/A3 comparators; improvements against those are209.525383% and419.429336%, respectively, not claims over a complete literature survey.
