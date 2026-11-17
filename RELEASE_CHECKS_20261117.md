# Release verification,17November2026

Passed before publication:

- All six85D–90D packing primary replays and all three independent audit batches, in the final portable repository package.
- Complete36D–48D scalar rank-two suite: all five compact angular gates, exact infinite-tail checks, every dimension-specific predecessor radial/profile dependency, exact dimension restriction and origin transfer, and final conservative rounding.
- Strongest46D quadratic Fourier/profile replay and separate exact degree123 Laguerre inversion.
- Strongest48D combined Fourier sum-of-squares/profile replay, separate degree83 Laguerre inversion, symbolic SOS reconstruction, origin PSD and exact rounding.
- The actual compressed-input wrapper materialized every file and checked both stored and decoded SHA256 values. All decoded proof inputs are byte-identical to their verified sources.
- Registered scalar values agree exactly between the package and root JSON registry. Python syntax checks and git whitespace/diff checks passed.

The independent Fourier/profile replays use1024-bit Arb. The scalar compact gates use60–100 decimal interval digits, with exact polynomial reconstruction and Sturm sign checks. A shared arithmetic backend is not an independently formalized theorem prover. Public source-minimum and positive-definite/Fourier theorems remain explicit mathematical inputs.

The unchanged unrestricted upper-bound suite and unrelated historical constructions were not rerun for this release. No GitHub Actions workflow is installed or dispatched. No harmonic-dependent matrix candidate or unproved93D coset claim is included.
