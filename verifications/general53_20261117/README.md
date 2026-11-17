# Certified unrestricted53D upper bound

Sphere-packing density is **<7.40755520e-8** in dimension53. This uses a degree151 positive-factor Fourier witness and universal annular neighbor-count budgets. The kernel is zero; no lattice assumptions enter.

Run `python run_all.py` for full1024-bit Fourier/tail/radial replay and independently authored exact cap/normalization checks. `--quick` checks exact algebra against the saved interval receipts, not a fresh sign traversal. Install `python-flint` and `sympy`; do not disable assertions. The hash-checked wrapper runs in a temporary copy.

[Independent theorem and audit scope](AUDIT.md). The full audit covered1309 annuli and1,810,331 nodes; shared FLINT/Arb is explicitly not a separate interval backend. The archived factor/root and annular certificates are preserved. Status strings in old discovery files are not proof claims.

The external comparison is the Afkhami-Jeddi, Cohn, Hartman, de Laat and Tajdini(2020) numerical spinless-bootstrap benchmark, approximately7.413767463149593e-8, not a certified upper bound or guaranteed-accuracy value. See the repository numerical benchmark table and its source caveats. This is not an exhaustive global-priority claim.
