# Independent 44D audit

PASS: Delta_44 < 0.000001842868599. This is a computer-assisted unrestricted-density certificate, not an optimality or worldwide-priority claim. See AUDIT.md for the mathematical review and precise independence scope.

## Replay

Python 3.11+ and python-flint 0.9.0 are sufficient. Run without optimization:

```sh
python verify.py
```

The entrypoint checks the exact frozen input hashes, replays all global source proof obligations, and executes separately implemented Fourier, angular PSD, origin Gram, radial coupling, and exact rounding checks. Its final result is `audit_result.json`; fresh logs are `complete_replay.log` and `independent_radial.log`. There is no optimization, network access, external data lookup, or repository-specific absolute path.

The complete source verifier already executes both its angular and cyclic-domain verifiers. The isolated independent script can also be run as `python independent_radial.py [package_directory]`.

The primary numerical inputs are frozen under `frozen_package/`. `INPUT_SHA256.json` covers all executable proof dependencies and witness files. The displayed comparator is the pinned previous upper 0.0000019017703144, not a lower bound. Lower-bound gap calculations belong to the separately sourced benchmark report and are not needed for this certificate.
