# Dimension five

Current strongest repository bound: **<0.510522525**. See the [radial-degree-four refinement](refinement_p4/), its [proof](refinement_p4/PROOF.md), and [independent audit](refinement_p4/INDEPENDENT_AUDIT.md). It removes approximately 3.71% of the updated paper-reported gap. Run `python refinement_p4/verify.py`.

The [earlier refinement](refinement_pd2/) proves **<0.510584101** and remains available.

## Preserved initial certificate

Frozen unrestricted sphere-packing density upper bound **0.511534153**, independently audited within this research effort. [Proof](PROOF.md) · [Audit](INDEPENDENT_AUDIT.md).

Run `python verify.py` with Python3.11+, python-flint0.9.0 and NumPy. Never disable assertions. No solver or network is needed for the new witness replay.

The proof uses four rigorous published local-neighbor-count bounds from arXiv:1206.2608 Table2 as external mathematical inputs; their original SDP certificates are not reproduced here. The audit consists of source review, independent reruns and objective reconstruction, not a clean-room implementation or proof-assistant formalization. No worldwide-priority or optimality claim is made.
