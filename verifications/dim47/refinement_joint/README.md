# Dimension 47: joint angular and radial refit

The target strict upper is **Delta_47 < 0.0000006010083731**. Check full_replay_result.json for the complete primary replay status and the separate independent audit for its status.

Within the repository run python verify.py for registration, hashes and full primary/independent chains. Standalone copies run python verify_primary.py and python independent/verify_all.py with Python3.11+ and python-flint0.9.0. No optimizer or network is required. Do not use Python -O.

The new target-dimensional degree-19 angular polynomial is unchanged from the preceding border certificate. The degree-111 radial polynomial is reoptimized for its improved angular cost, using a degree-31 Laguerre correction and exact root-guided sign repairs. Every Fourier, global sign, pair, PSD, cyclic and rounding gate is replayed.

This is a computer-assisted mathematical certificate, not a proof-assistant formalization, external peer review, optimality result or exhaustive priority claim. Full discovery/export provenance is under discovery/.
