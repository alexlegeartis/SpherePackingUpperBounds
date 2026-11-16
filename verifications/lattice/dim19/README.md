# 19D lattice-only density bound

**Density ≤ 0.007808430327.** This is an unconditional lattice bound, not an unrestricted sphere-packing bound.

The certificate combines Zubrilina's modulo-two class framework with global short-class angular coupling. It requires no prior density bound, trace gap, determinant lemma or arithmetic exclusion. [METHOD.md](METHOD.md) explains the radial/angular combination; [ANGULAR_PROOF.md](ANGULAR_PROOF.md) proves the continuous angular counting and exact projective caps. The certificate uses the quadratic and quartic caps only, not the later exploratory generalization in the proof document.

Run `python run_all.py` with Python 3.10+ and SymPy. `python-flint` is an optional public SymPy acceleration dependency; without it the same exact checks remain correct but may take substantially longer. No private interpreter path is required. Do not use `python -O`. Both the original exact verifier and a separately written independent verifier run. They reconstruct the Fourier transform, check all global signs, all 120 radial intervals, all three continuous angular cuts, every class budget, and the exact final normalization. No optimizer, network or research-directory import is needed. The independent verifier writes a local JSON audit report. The compressed release runner executes in a temporary directory.

Certificate SHA256: `7196d4c4c1b0e33896413ca2ef064cef6c0ffecd77049b49532010964ec07447`.

## Comparison and attribution

The earlier named upper bound is 0.00782447361 from Nina Zubrilina, *Improved lattice-packing bounds in small dimensions*, https://arxiv.org/abs/2610.35830 . `comparison_certificate.json` is the original `certificates/n019.json` from https://github.com/AstroNinja389/lattice-packing-bound/tree/340a0e93bb95942c2b1b80092be43149504970eb . It is preserved only as an attributed benchmark, not an assumption of the new proof. Its source verifier is available in `../prior/`; optional replay: `python ../prior/taylor.py comparison_certificate.json`.

This is an improvement over that named comparison, not an optimality or exhaustive worldwide-priority claim. Delsarte/Gegenbauer positivity and projective relative bounds are established tools; the polynomial-bound family is not claimed novel here.
