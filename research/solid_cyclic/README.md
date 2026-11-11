# Solid-harmonic cyclic search

The final 17D and 18D cyclic witnesses have been independently audited and registered in ../../verifications/dim17/ and ../../verifications/dim18/refinement_cyclic/. Earlier rejected and numerical-control witnesses are retained for reproducibility and are not certified results. In particular, 16D, 22D and 23D numerical outputs here are not new registered bounds.

Install requirements.txt for discovery. The root verification dependencies are sufficient for exact replay. Run:

    python search.py --dimension 17 --degree 31 --rounds 5 --output candidate.json

Search parameters and round histories are retained in candidate17_v2.json. The exact final witness certificate17.json, rather than a solver status, defines the proof. This research snapshot preserves its original status metadata; the registered package reconciles it with the completed independent audit.

The far-domain failure of the first candidate is preserved in rejected17_before_full_separation.json. It motivated adding all real angular stationary points to numerical separation. Full compact, far and pair interval checks remain mandatory. The cyclic method and solid-harmonic factors are explained in PROOF17.md; it must not be described as satisfying the original pointwise-sign Theorem 1.1 without proof.

For candidate export, use the dimension-generic rationalizer: `python ../annular_three_point/rationalize.py candidate.json witness`. Then test with `python verify.py witness/dim17.json`. Export is unverified; the wrapper may reject it, and numerical reruns can vary. The retained registered exact certificates are the reproducible proof inputs.
