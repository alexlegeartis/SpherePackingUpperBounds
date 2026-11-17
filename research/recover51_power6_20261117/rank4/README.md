# Four-profile compact-kernel search in 51D

**Numerical research only. No bound or record is established here.**

The first family is `[4,8,12,16]`, enlarging the previously published `[8,12,16]` feature space. Fixing the original exact Fourier-positive degree-71 scalar keeps certification of its Fourier sign available while testing whether additional local-kernel freedom materially reduces the origin cost. A strictly positive numerical PSD margin means exact inclusion of the lower-rank face is approximated, not imposed.

There are thirty 4-by-4 harmonic blocks and a 5-by-5 origin block. Optimization uses the shifted Legendre basis through degree three, then exports monomial-profile matrices. For `u=beta^spacing`, the ten support-crossing controls use the exact cubic Bernstein transform

`T[j,i] = binomial(j,i)/binomial(3,i)` for `i<=j`, zero otherwise.

The identity `B_3(u) T = (1,u,u^2,u^3)` was checked symbolically. Exact rational angular derivative-root isolation supplies critical-point cuts. Triangle separation includes initial boundary grids, all twelve edges with dense sampling and one-dimensional local extrema, all six faces with two-dimensional differential evolution, and three-dimensional differential evolution plus random tests. Those are discovery procedures, not continuum sign proofs.

Source: published 51D rank-three power-eight certificate at commit `34f59648d034e92a52768c95af96221a8f589ac9`, archive SHA-256 `7d1b48627948429301d916f3ed04beb255d33c7cd325cb47a7ca81acbb7128e6`. Adjacent `../bootstrap.py` reconstructs `../base` from that archive.

Run from a checkout with that source and numpy, scipy, sympy, python-flint, cvxpy:

```sh
OPENBLAS_NUM_THREADS=1 python search.py --power 4 --spacing 4
```

Each iteration saves the full numerical candidate and cut lists in a separate family subdirectory. `unscreened*.json` is saved before separation; `candidate*.json` includes screen outcomes. A successful screen still needs exact rational repair, complete triangle-domain and radial coverage, all tail/counting inputs, and independent replay.

## First numerical outcome

The first `[4,8,12,16]` solve returned origin cost about 1.34825, worse than the published 1.32396, and failed both edge and angular critical-point screens. This run was stopped. A possible conditioning issue is that adding a new leading feature creates zero baseline Bernstein controls; imposing a uniform strict margin then excludes the baseline face substantially.

The next tested family is `[8,12,16,20]`. It preserves the original leading feature, and the old quadratic Bernstein controls become strict cubic controls by ordinary degree elevation. This separates additional feature freedom from the vanishing-leading-control effect. No successful result is asserted yet.
