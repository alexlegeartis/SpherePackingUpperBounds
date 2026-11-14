"""Replay the aggregate's final 48D rounding and shared origin checks separately.

The initial historical aggregate stops before these final checks because one
pre-existing 16D registry assertion failed. This records their actual replay;
it does not relabel the initial aggregate as successful.
"""
from fractions import Fraction
from pathlib import Path
import json,subprocess,sys

if not __debug__:raise RuntimeError('Assertions must remain enabled')
ROOT=Path(__file__).resolve().parents[2]
a=json.loads((ROOT/'verifications/dim48/continuum_three_point/independent_cyclic40_radial.json').read_text())
assert Fraction(a['density_upper'])<Fraction('4.13385739314e-7')
for p in ['verifications/origin_slack/verify_algebra.py','verifications/cyclic_origin_slack/verify_algebra.py']:
    subprocess.run([sys.executable,str(ROOT/p)],cwd=ROOT,check=True)
print('PASS final 48D rounding and both shared origin algebra checks')
