# Dimension four

Current strongest repository bound: **<0.633430327**. See the [radial-degree-six refinement](refinement_p6/), [proof](refinement_p6/PROOF.md), and [independent audit](refinement_p6/INDEPENDENT_AUDIT.md). It closes approximately 13.88% of the updated paper-reported gap. Run `python refinement_p6/verify.py`. Earlier certificates below remain preserved.

## Earlier radial-degree-four refinement

Earlier bound: **<0.634526023**. See the [radial-degree-four refinement](refinement_p4/), its [proof](refinement_p4/PROOF.md), and [independent audit](refinement_p4/INDEPENDENT_AUDIT.md). It closes approximately 8.19% of the updated paper-reported baseline gap. Run `python refinement_p4/verify.py`.

## Preserved earlier certificate

Independently audited unrestricted density upper bound **0.635766243**. See [proof and external assumptions](PROOF.md) and [audit](INDEPENDENT_AUDIT.md). This closes approximately 1.77% of the pinned upper–lower gap, not 20%.

Run `python verify.py` from this directory. Python 3.11+, python-flint 0.9.0 and NumPy are required. No solver or network is needed. Do not disable assertions. The runner uses a temporary directory and leaves the frozen witness unchanged.

The four published near-neighbor count bounds are external mathematical inputs. The new witness and its global inequalities are included in full. Independent computational audit is not external peer review, formal proof, or a priority claim.
