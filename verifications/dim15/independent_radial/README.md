# Independent 15D curved-profile radial audit

PASS for every assigned radial, normalization, origin-PSD and exact rounding
gate. The resulting density expression is strictly below **0.033671603**.
The sphere-packing conclusion additionally uses the matching full independent
angular/cyclic proof identified in `ANGULAR_DEPENDENCY.json`.

Run `python verify_all.py` with Python 3.11+ and python-flint 0.9.0. Do not use
Python -O. No optimizer, network, primary verifier or sibling folder is needed
to replay the radial audit itself.

`AUDIT.md` contains the proof review, exact odd-dimensional volume formula,
normalization, full curved pair argument and theorem linkage. `INPUT_SHA256.json`
pins all mathematical inputs. `SHA256SUMS` covers the portable package.
`audit_result.json` and `full_replay.log` record the successful replay.
