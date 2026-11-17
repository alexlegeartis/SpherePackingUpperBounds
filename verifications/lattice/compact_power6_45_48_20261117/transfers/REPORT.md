# Independent45–47D restriction/scalar audit

**COMPLETE PASS: the source48D main geometry has now passed full independent 384-bit replay, with 183,451 nodes and 91,726 leaves.**

The only source kernel has SHA256169ea8d0e2ec662e7d450641894d1e3ba5f280623d2d6a659d1308401d45a728; see exact_result JSON files for the authoritative full SHA and exact frozen target hashes. No earlier rejected kernel is used.

The checker constructs each source angular matrix polynomial directly in the48D normalized Gegenbauer basis, reconstructs the corresponding target polynomial from the frozen target matrices, and verifies exact coefficientwise equality to a positive scalar multiple of the source. All target matrices and origin Schur complements are checked exactly. This proves transfer of the entire plane sign theorem using the now-completed source theorem; it is not merely a numerical restriction or a claim based on scalar Gegenbauer coefficients alone.

Each target additionally has its own degree71 positive-factor Fourier polynomial, exact Laguerre inversion,633 contiguous annuli with independently recomputed mod2 caps and exact budgets, and71-root negative-tail verification. All physical coefficients and kernel data are linked exactly to the frozen scalar inputs. Fresh isolated radial and tail replays have been run, using the reviewed producer Arb/Taylor implementation rather than falsely claiming a new interval backend.

Verified conservative upper bounds:45D1.1728551765e−6,46D7.721596228e−7,47D5.130033534e−7. Even46D normalization is rational. Odd45D/47D normalization is squared and verified using a rigorously directed rational Machin lower bound for pi. All dependencies of these three bounds have now passed the described checks. The shared source trace SHA is `2145900b7ea753d5d711fd896ebca614961a56ae9427bf8b6fe1623185c3d156`.
