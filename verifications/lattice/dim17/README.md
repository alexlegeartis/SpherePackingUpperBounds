# 17D lattice-only density bound

**Density ≤ 0.013832417384.** This is an unconditional lattice bound, not an unrestricted sphere-packing bound.

The certificate combines Zubrilina's modulo-two class framework with global short-class angular coupling. It requires no prior density bound, trace gap, determinant lemma or arithmetic exclusion. [METHOD.md](METHOD.md) explains the radial/angular combination; [ANGULAR_PROOF.md](ANGULAR_PROOF.md) proves the continuous angular counting and exact projective caps. The certificate uses the quadratic and quartic caps only, not the later exploratory generalization in the proof document.

Run `python run_all.py` with Python 3.10+ and SymPy. `python-flint` is an optional public SymPy acceleration dependency; without it the same exact checks remain correct but may take substantially longer. No private interpreter path is required. Do not use `python -O`. Both the original exact verifier and a separately written independent verifier run. They reconstruct the Fourier transform, check all global signs, all 120 radial intervals, all three continuous angular cuts, every class budget, and the exact final normalization. No optimizer, network or research-directory import is needed. The independent verifier writes a local JSON audit report. The compressed release runner executes in a temporary directory.

Certificate SHA256: `4eeb6ebef7e74c1f4a5d26af00f949411cfb57a8a468516cdf644778d346b626`.

## Comparison and attribution

The earlier named upper bound is 0.0138839313 from Nina Zubrilina, *Improved lattice-packing bounds in small dimensions*, https://arxiv.org/abs/2610.35830 . `comparison_certificate.json` is the original `certificates/n017.json` from https://github.com/AstroNinja389/lattice-packing-bound/tree/340a0e93bb95942c2b1b80092be43149504970eb . It is preserved only as an attributed benchmark, not an assumption of the new proof. Its source verifier is available in `../prior/`; optional replay: `python ../prior/taylor.py comparison_certificate.json`.

This is an improvement over that named comparison, not an optimality or exhaustive worldwide-priority claim. Delsarte/Gegenbauer positivity and projective relative bounds are established tools; the polynomial-bound family is not claimed novel here.
