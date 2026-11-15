"""Lemma 969 (cyclicslack, D4-12): the accepted raw V7 kernel K by polarization
of the retained evaluator, the physical cyclic term Gamma, and the allocated
slack H_R at an anchor. Run from any folder; raw_v7.py sits beside this file.

The anchored energy is E(C) = A + 2 sum_x h(x) + sum_{x,y in C} K(x,y), a
quadratic form in the star (audit section 1). raw_v7.evaluate returns
outward rational enclosures of E (complete: every component and the gauge).
Polarization: for a != b,
    K(a,b) = (E({a,b}) - E({a}) - E({b}) + A) / 2,
which depends only on (|a|^2, |b|^2, 2 a.b). K = 0 if |a| or |b| > 7/4.
Cyclic term of the triangle (0,y,z): Gamma = K(y,z) + K(-y,z-y) + K(-z,y-z).
H_R at anchor 0 = sum over ordered pairs y != z of the star with |y|,|z| <= 7/4
of -Gamma(y,z)/k_R, k_R = 3 if |y-z| <= 7/4, else 1 (plans D4-12).
"""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import importlib.util, sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("raw_v7", HERE / "raw_v7.py")
raw = importlib.util.module_from_spec(spec); spec.loader.exec_module(raw)
A = raw.A
R2 = Q(49, 16)
O = (Q(0),) * 4
dot = lambda x, y: sum(a * b for a, b in zip(x, y))
sub = lambda x, y: tuple(a - b for a, b in zip(x, y))

def E(vectors):
    return raw.evaluate([O] + list(vectors), O)["total"]

@lru_cache(None)
def E1(m, vec):
    return E([vec])

_K = {}
def K(a, b):
    """Enclosure (lo, hi) of K(a,b), a != b, both nonzero, |a-b| >= 1."""
    m, n, h = dot(a, a), dot(b, b), 2 * dot(a, b)
    if m > R2 or n > R2:
        return (Q(0), Q(0))
    key = (min(m, n), max(m, n), h)
    if key not in _K:
        e2 = E([a, b]); ea = E1(m, a); eb = E1(n, b)
        _K[key] = ((e2[0] - ea[1] - eb[1] + A) / 2, (e2[1] - ea[0] - eb[0] + A) / 2)
    return _K[key]

def Gamma(y, z):
    t = [K(y, z), K(tuple(-c for c in y), sub(z, y)), K(tuple(-c for c in z), sub(y, z))]
    return (sum(x[0] for x in t), sum(x[1] for x in t))

def H_R(star):
    """Enclosure of H_R at anchor 0; star = vectors with |v|<=7/4. Returns (lo, hi, per-key breakdown)."""
    lo = hi = Q(0)
    for i, y in enumerate(star):
        for z in star[i + 1:]:
            g = Gamma(y, z)
            assert g[0] <= 0 or g[1] <= 0
            k = 3 if dot(sub(y, z), sub(y, z)) <= R2 else 1
            lo += 2 * (-g[1]) / k; hi += 2 * (-g[0]) / k
    return lo, hi
