# Global-positivity follow-up

The initial Fourier sample-sign obstruction has been repaired by imposing Bernstein coefficient signs on a partition of the finite interval and Taylor coefficient signs on the complete infinite tail. Numerical arithmetic remains discovery only; exact Sturm verification independently certifies the rounded rational polynomials.

`search_global.py` additionally imposes angular Bernstein coefficients at corners of every fine radial rectangle, extending each angular interval to the upper-neighbor cosine. This matches the continuum interpolation proof used by `verify_global.py`. These conditions cover full angular intervals, not only sampled angles. The annular coupling is sampled in discovery and verified separately over the full interval using Arb.

## Completed rigorous certificates

For the degree-23 run `candidates_global_full23.json`:

- Dimension 17: exact witness `witnesses_global_full23/dim17.json`, density <0.018783841.
- Dimension 18: exact witness `witnesses_global_full23/dim18.json`, density <0.013768221.

Both pass exact PSD, exact Fourier identity, exact Fourier positivity and exterior negativity, all 153 angular corner intervals, and the complete annular inequality. The exact witnesses and objective expressions define the result; rounded candidate densities are not the proof.

These are **not benchmark improvements**: the pinned Cohn table lists 0.0184640903350649 and 0.0134853404450862 respectively. They should remain research/regression certificates rather than being promoted to the main result registry. The index snapshot is not a worldwide-priority audit.

The 22D/23D degree-23 full-Bernstein run failed exact PSD and is rejected. Its solver objective must not be reported as a bound.

## Proof specialization

Use Theorem 1.1 of Cohn–de Laat–Salmon (https://arxiv.org/abs/2206.15373) with the rational exclusion/truncation squares recorded in the witness and five radial hat nodes. Define the normalized Fourier-eigenbasis functions using L_k^(n/2-1)(2X)/binomial(k+n/2-1,k). The exact verifier reconstructs the Fourier identity in rational arithmetic, with n/2 kept rational also for odd dimensions.

The kernel construction, isolated origin feature, PSD addition theorem, angular-corner interpolation, and annular Taylor enclosure are exactly those described in `reference/PROOF.md`, with dimension, degrees, radii and matrices replaced by this witness. The interval cosine monotonicity hypothesis R²<2r² is checked. The density objective is

    pi^(n/2) / Gamma(n/2+1) * (r²/4)^(n/2) * (P(0)+A0[0,0])/H(0).

It is enclosed using 256-bit Arb arithmetic. Both verification scripts are required; one alone is not the complete certificate. The underlying continuum theorem and positive-definite Gegenbauer addition theorem are mathematical inputs, not formalized here.

## Completed fine degree-31 run

Both fine-grid witnesses pass all local computational gates. Dimension 17 yields <0.018582096, still above the pinned comparator. Dimension 18 yields **<0.013408931**, about 0.5666% below the pinned comparator. See `PROOF_DIM18.md`, `witnesses_global_fine31/dim18.json`, and `audit_global_fine31_18.log`. Independent review passed; see ../../verifications/dim18/REVIEW.md. Updated searches support the pinned comparison but do not establish worldwide priority. The discovery run has finished; no solver remains active.

## Reproduce

Install requirements.txt. From this directory, run `python search_global_fine.py --dimensions 17 18 --degree 31 --output candidates_global_fine31.json`, followed by `python rationalize.py candidates_global_fine31.json witnesses_global_fine31`. Discovery outputs remain unverified and may vary by solver/platform. The retained exact witnesses are the reproducible inputs. Run `python run_certificate.py witnesses_global_fine31/dim18.json` to replay both producer checks; the registered independent path is `../../verifications/dim18/verify.py`. The rationalizer chooses its own unverified objective threshold; the retained final18D witness uses the tighter independently checked strict decimal 0.013408931.
