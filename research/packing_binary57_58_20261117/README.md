# Verified binary complex-multiplication extensions in dimensions 57 and 58

Run `python run_all.py` for exact setup, parent-source/basis replay, completed-log integrity and 40 independent brute-force controls. Run `python run_all.py --full` to additionally compile and replay BOTH exhaustive exclusions (roughly 25 minutes, hardware dependent). Requirements: Python 3, numpy, python-flint and g++ with IEEE binary64 and long double with at least 64 significand bits. Do not use Python -O or unsafe floating-point compiler flags.

Copied source inputs and proof scripts remain intact. The wrapper works in a temporary copy so generated reports do not modify frozen inputs. Quick mode does not reprove exclusion by computation: it checks the archived completed logs. The primary setup script retains its historical independently-pending wording; final independent status is PASS in cm57_order2_independent_20261117/AUDIT.md.

See cm57_deeper_20261117/BINARY_CM_EXTENSION.md for construction and cm57_order2_independent_20261117/AUDIT.md for the independent proof. The unchanged parent minimum, source-section and basis provenance are inherited from cm57_independent_20261117, whose complete exact input is included. The original denominator-three files are dependencies only, not the principal external comparison. No exploratory screens or optional local-extremality claims are included.

Rigorous ordinary density floors: 57D >= 1.61599884109140e-11; 58D >= 5.58243739609337e-12.

External comparisons are derived source consequences, not printed dimension-specific records: 57D MW56 times A1, Schuett-Shioda Theorem 13.35 at p=29, squared center density (5/2)^57/(10*29^8), gain 82.574185835055%; 58D Gamma72 perpendicular to a rank-14 Leech block, Nebe 2012 section 4.2 plus exact section verification, squared center density 2^36/3, gain 24.749001701336%. No exhaustive global-priority claim is made.

## External comparisons

The57D comparison is the minimum-matched MW56 × A1 construction from Schütt–Shioda, Theorem13.35(p=29), with squared center density(5/2)^57/(10*29^8). The58D comparison is the Γ72 section perpendicular to a primitive rank14 Leech block, derived from Nebe(2012), Section4.2; the existing repository's [exact section verifier](../novel_routes_20261113/audits/verify_gamma72_leech_sections.py) gives squared center density2^36/3. These are consequences of published inputs, not dimension-specific numerical records printed by those sources. [Nebe paper](https://arxiv.org/abs/1008.2862) · [Schütt–Shioda source](https://web.archive.org/web/20200209144320id_/http://www2.iag.uni-hannover.de/%7Eschuett/BookMWL17.pdf).

Exact squared-density ratios give increases82.57418583505537% and24.74900170133646%, respectively. Our own order-three and previous half-coset constructions are not the principal comparators. RESULTS.json records exact determinants, center-density squares, comparison squares and rigorous occupied-density floors. These are the strongest attributable comparators retained in the documented source review, not an exhaustive global-priority guarantee.
