# Origin-coordinate slack research

Exact Schur-complement postprocessing of the existing 4D, 5D, 19D, 20D and 21D witnesses. This is ordinary PSD slack removal, not a new multipoint theorem or a worldwide-priority claim.

**Research only: RESULTS.json is unchanged.** The original continuum certificates are inherited and have not been replayed by these scripts. All derived density claims remain conditional on those source checks and external mathematical inputs.

Source: alexlegeartis/SpherePackingUpperBounds commit 85e395353a937369e5cff6400f66337967812270. Files dim4_witness.json and dim4_psd.json correspond to the p6 refinement; witness.json and witness_psd_exact.json to the 5D p4 refinement. The other source witnesses correspond to the respective dimension packages.

Install requirements.txt, then run:

    python origin_slack.py
    python independent_verify_5d.py

The first command checks exact SymPy LDL and residual factorization in five dimensions. It rewrites origin_slack_results.json. The second independently checks the 5D factorization and reconstructs its tightened objective using Fraction and Arb. Do not use Python -O. See ORIGIN_SCHUR_PROOF.md for the unchanged-continuum-constraints argument and normalization.
