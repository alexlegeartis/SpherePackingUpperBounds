# Certified unrestricted upper bounds in51D and52D

**51D density<1.52983244e-7.52D density<1.06516626e-7.** These are unrestricted sphere-packing results, not lattice-only bounds.

The reusable method is a positive Fourier witness plus a universal per-center neighbor budget.51D retains the previously audited cyclic PSD angular kernel;52D sets the kernel to zero. A sharp exact spherical-code polynomial controls the first thin annulus, while exact Levenshtein-type polynomials control subsequent annuli. No coset, lattice or antipodal premise enters these neighbor caps.

## Replay

Run `python -B run_all.py` with numpy, mpmath and python-flint installed, without `-O`. The hash-checked wrapper uses a temporary copy and leaves the published package unchanged. `--quick` verifies exact Fourier/cap/PSD/density algebra using the recorded interval receipts; it does not freshly recompute interval signs. Full mode replays all global Fourier/tail and radial gates, plus51D's angular gate.

The independent full1024-bit runs passed51D633 annuli/1,099,429 nodes and52D909 annuli/1,706,071 nodes. Separately authored checks reconstruct Fourier transforms by Laguerre recurrence, prove contact-polynomial signs by exact Sturm sequences, reconstruct later cap coefficients, check interval coverage and recompute safe bounds with rational arithmetic. Shared FLINT/Arb interval replay is not a separate arithmetic backend. All inputs and complete logs are included; see THEOREM_REVIEW.md for the theorem, normalization and assumptions.

## Attributable numerical comparison

Afkhami-Jeddi,Cohn,Hartman,de Laat,Tajdini(2020), numerical spinless bootstrap, supplies the external comparison values. The highest supplied orders areN457 in51D andN466 in52D, giving approximate density benchmarks1.5321362440043e-7 and1.0661233909836e-7. Our safe decimals are0.1503654791% and0.0897767549% smaller, respectively. These are modest improvements relative to numerical estimates, not reductions of certified prior upper bounds. The source spectra are not globally sign-certified Cohn–Elkies functions, and the same-radius lift is conjectural in general.

[Complete49–96 numerical benchmarks and source provenance](../../NUMERICAL_BENCHMARKS49_96.md). No global optimality or exhaustive priority claim is made. No new53D result is included, and the lattice-only table remains restricted to dimensions at most48.
