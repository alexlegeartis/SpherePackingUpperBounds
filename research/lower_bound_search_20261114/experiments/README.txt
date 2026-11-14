OPTIONAL NUMERICAL EXPERIMENTS

These scripts perform discovery, not rigorous global minimum verification.
Historical outputs under ../results/ are kept separately from regenerated
files in this directory. Install experiments/requirements.txt first.

The exact finite proofs and the exact cyclic-code 34D local-optimality
certificate need only NumPy and live one directory above this folder.

The separate minimum-twelve self-orthogonal-fiber search and its exact audit
are documented in ../NORMALIZATION12.txt. Its driver selforth12.py requires
NumPy and a C++ compiler for search; auditing saved generators needs only NumPy.

Examples, run from the research package root:

    OPENBLAS_NUM_THREADS=1 python experiments/replay.py KP34
    OPENBLAS_NUM_THREADS=1 python experiments/replay.py CTturyn1_scale_vertex
    OPENBLAS_NUM_THREADS=1 python experiments/replay.py from38_section37_0
    OPENBLAS_NUM_THREADS=1 python experiments/replay.py new36

The last command regenerates the new36 contact file used by:

    OPENBLAS_NUM_THREADS=1 python experiments/section_metric_screen.py

The alternative CT vertex's ten exact contact-preserving directions can be
regenerated after its replay command:

    OPENBLAS_NUM_THREADS=1 python experiments/ct1_vertex_flex.py

For its ten-parameter vertex screen, first run these from this directory:

    ./enum_shell 36 10.0000001 1500000000 CTturyn1_scale_vertex_R.txt \
        CTturyn1_vertex_short10.txt CTturyn1_vertex_short10.log
    OPENBLAS_NUM_THREADS=1 python ct1_vertex_screen.py

The supplied CTturyn1_vertex_flex.npz stores integral Q at denominator 16,
ten integer directions D, and checked contact rows. The scale vertex NPZ's Q
is the physical rational Gram, represented exactly in binary at denominator
16. Other integer seeds use the normalizations in the main README and logs.

The tensor/code and section experiments do not require any remote services.
The basis and minimum of the original catalogue seeds are literature inputs
unless explicitly verified by one of the exact scripts above this directory.

Some experiment scripts overwrite their generated files, including numerical
candidates and regenerated kernel arrays. Use a copied experiment directory
when preserving an existing run matters. A BKZ basis minimum is only an upper
bound on the true lattice minimum. Floating enumeration completion supplies
neither a rigorous minimum certificate nor a record claim.
