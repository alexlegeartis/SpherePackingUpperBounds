# Independent 18D review, 10 November 2026

Both producer full-domain replay and independent verification pass for the frozen certificate with SHA256 286d466bfd83e3ec8543f7fa7312655432701289a51f84c76a369afbdf0ee137.

The independent implementation is adapted from the repository's independent 19D verifier, not from the new candidate verifier. It uses rational Schur elimination, the direct monomial Fourier transform formula, explicit normalized Gegenbauer coefficients, integer-square-root cosine caps, independent Machin bounds for pi, and independently reconstructed Arb Taylor enclosures. It checks 17 PSD blocks with 86 positive pivots, 561 angular polynomials, and 34 annular leaves. It certifies density strictly below 0.013408931.

The exact exclusion squared is 9/4. The exact truncation squared is 26999999999999997/10000000000000000, not 27/10. Both scripts use the stored rational. The former may be written r=3/2; the latter must not be silently rounded in a proof.

Theorem 1.1 of Cohn, de Laat and Salmon, https://arxiv.org/html/2206.15373v1, was read directly. Set f2=P(pi|x|^2)exp(-pi|x|^2)/H0 and divide the entire radial-hat kernel by H0. The Fourier transform has polynomial H, with H(0)=H0>0. The scalar function is Schwartz. Kernel PSD follows from rational PSD matrices and the normalized Gegenbauer addition theorem. The isolated origin coordinate and closed annular support are permissible: the theorem does not demand kernel continuity.

The angular inequality is global, not sampled. On each fine radial rectangle the kernel is bilinear in squared radii for fixed cosine, hence a convex combination of four corner polynomials. The allowed cosine (u+v-r2)/(2 sqrt(uv)) is nondecreasing in each radius because R2<2*r2. The certified cap at each corner's upper-neighbor radii therefore covers all cells using that corner. Every required cosine lies in the verified interval [-1,cap]. The pair inequality is verified throughout each coarse radial interval by exhaustive Taylor subdivision. The physical tail and Fourier positivity hold on their complete half-lines by exact Sturm counts. No neighbor-count theorem or lattice-only assumption is used.

Normalization is independently checked: the density multiplier in dimension 18 is pi^9/9! times (r2/4)^9, and objective (P(0)+A0[0,0])/H0. The final upper is proved using a rational upper bound for pi, not merely a decimal printout.

Literature comparison: Cohn's indexed general-packing table gives 0.0134853404450862, attributed to Afkhami-Jeddi et al., https://arxiv.org/abs/2006.02560. The newer https://arxiv.org/html/2512.04829v2 Table 1 covers 4–7 and 9–16, not 18. Refreshed searches found no stronger finite 18D unrestricted bound. Live Cohn retrieval fails and the available index was crawled four months earlier; this establishes improvement against that pinned comparator, not worldwide priority. Relative upper decrease using outward 0.013408931 is 0.5666111686 percent.

Replay: PYTHONPATH=/tmp/density_packages python independent_verify.py. The temporary PYTHONPATH is environment-specific; install python-flint normally elsewhere. Assertions must be enabled.
