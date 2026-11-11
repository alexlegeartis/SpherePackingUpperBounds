# Independent 46D extended-pair certificate audit

PASS: Delta_46 < 0.0000008777229084.

Run `python verify_all.py` using Python 3.11+ and `python-flint==0.9.0`.
Do not use Python -O. No optimizer or network is required.

This package replays the full angular/cyclic proof and independently checks the
new radial export, exact 100 ppm M guard, complete Fourier ray, physical ray
from X = 40, every point of the full coupled interval [6*pi, 40], PSD origin
extension and exact rational density rounding. It also proves P(31.018) > 0,
which places this witness outside the old forced-negative tail from 7*pi.

See AUDIT.md for the theorem, scope and proof review. INPUT_SHA256.json pins
all frozen inputs. SHA256SUMS covers the portable package. audit_result.json
and full_replay.log record the aggregate replay.
