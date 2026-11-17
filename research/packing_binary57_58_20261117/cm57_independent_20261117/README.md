# Independently verified MW56 coset and57D/58D lattices

See AUDIT.md and AUDIT_RESULT.json for completed scope. This package certifies coset minimum20/3, lattice minimum10, and density floors57D1.53307111021240e-11 and58D5.02419365648403e-12. No global-record claim is made.

Run check_setup.py with python-flint and NumPy, parent/verify_parent.py under input/, and check_density.py. Compile and run independent_search.cpp as documented in AUDIT.md. The independent complete enumeration took about16minutes and3.65b nodes; reaching its6b cap would fail, not certify absence. The control executable changes only the admissible dimension range for40 brute-force comparisons. The original producer and independent logs are both preserved.
