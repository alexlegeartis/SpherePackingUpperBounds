# Verified 81-dimensional binary-coset extension of a Bachoc-Nebe lattice

Minimum 8, determinant 2, exact center-density square 2^80, ordinary density at least 2.87777416803990e-17. The reconstructed 80D parent is a Bachoc-Nebe ten-E8-block neighbor; no identification with a particular catalogue basis is required.

Run `python run_all.py` for exact parent/source-hypothesis reconstruction, input maps, local finite-search coverage, contact, Gram/determinant and rigorous density checks. `python run_all.py --full` additionally recompiles and reruns BOTH complete integer-only searches plus 3 positive controls. Requirements: Python 3, numpy, python-flint, and g++. Assertions must remain enabled. Inputs/scripts are frozen and checked by SHA256; replays work on temporary copies.

The full primary and independently written searches each exhaust all 16 glues, 5,877,535,660 nodes and 100 residue leaves. The independent implementation changes nearest-point computation, variable/tie order and carry dynamic program. Its final status is PASS. The primary source retains its historical PENDING_AUDIT result label; this is resolved by the independent audit, not silently edited.

The parent minimum 8 theorem from Bachoc-Nebe (1998), Theorem 4.7, remains an external mathematical premise. Source and exact constructive hypotheses are documented in the parent PROOF.md. CASE101_PROOF.md explains the local shell coverage, global code join, norm congruence and extension theorem. No optional local-extremality claim is included.

Comparisons: +3.1054957484% over the pinned preceding repository result; +80.3894724750% over the derived MW80×A1 example. The former is internal progress; the latter is an attributable external construction, not a claim that it is the strongest possible consequence of previous literature.
