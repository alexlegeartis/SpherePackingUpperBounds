# Extended pair search and exact export

The new model keeps minimum squared distance 6, Gaussian rate17/8 and the same verified angular polynomial. It replaces the stronger physical constraint P(X)<=0 from 7pi with the actual pair inequality through X=40, followed by P(X)<=0 on[40,infinity).

The radial span is a nonnegative multiple of the normalized degree 111 source plus degree 47 normalized Laguerre corrections. The first sampled candidates failed Fourier signs or missed a narrow positive pair lobe nearX=31.018. Exact-root cuts and this additional pair point produced the frozen coefficients. Every rejected checkpoint is retained in history/; those rejection files do not describe the final witness.

The final export multiplies M by exactly 10001/10000, a 100ppm guard. The complete interval pair proof, not the sampled peak scan, validates the guard. No optimizer result or numerical peak search is used in the final density proof.

radial_search.py --dimension 46 --degree 47 reproduces the search design and saved cuts using NumPy, SciPy and mpmath. audit_radial.py radial_perturb46_degree47.json reproduces rational export and radial checks using python-flint. A new optimizer replay need not produce identical floating coefficients. Exact export has been replayed; no new optimizer replay is implied.
