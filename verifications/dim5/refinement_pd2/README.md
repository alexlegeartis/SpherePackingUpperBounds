# Dimension5

Frozen unrestricted sphere-packing density upper bound **0.510584101**, independently audited within this research effort. [Proof](PROOF.md) · [Audit](INDEPENDENT_AUDIT.md).

Run `python verify.py` with Python3.11+, python-flint0.9.0 and NumPy. Never disable assertions. No solver or network is needed for the new witness replay.

The proof uses four rigorous published local-neighbor-count bounds from arXiv:1206.2608 Table2 as external mathematical inputs; their original SDP certificates are not reproduced here. The audit consists of source review, independent reruns and objective reconstruction, not a clean-room implementation or proof-assistant formalization. No worldwide-priority or optimality claim is made.
