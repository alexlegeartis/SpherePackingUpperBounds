"""Rigorous outward fixed-binary tensor Bernstein conversion.

The input is an integer power-coefficient tensor divided by a positive
integer denominator. Bounds must be exact, with 0 <= left < right on every
axis. The result encloses EVERY Bernstein coefficient, not sample values.
No sign decision, subdivision rule, or certificate status is implemented.

All arithmetic before the explicit int64 handoff uses Python integers and
FLINT integer matrix products. No floating-point number is accepted.
"""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import hashlib
import importlib.util
import time

import numpy as np


REFERENCE_SHA256 = '3e7a56d77a9f481b21ed65ab2d580c0a40be9c52745440abf7a958189e989631'


def _reference_path():
    # Active workspace and portable research-archive layouts are allowed.
    # The first existing candidate must have the exact pinned bytes. A stale
    # or modified active dependency cannot be silently bypassed by fallback.
    here = Path(__file__).resolve()
    candidates = [here.parent.parent/'base/new_work/pointed/exact_bernstein_fast.py']
    candidates.extend(a/'historical_snapshot/new_work/pointed/exact_bernstein_fast.py'
                      for a in here.parents)
    seen = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        if candidate.exists():
            if not candidate.is_file():
                raise ValueError('Exact transform dependency is not a file.')
            if hashlib.sha256(candidate.read_bytes()).hexdigest() != REFERENCE_SHA256:
                raise ValueError('Exact transform dependency hash mismatch.')
            return candidate
    raise FileNotFoundError('Pinned exact Bernstein transform dependency not found.')


REFERENCE_PATH = _reference_path()
_spec = importlib.util.spec_from_file_location('_interval_bernstein_exact_reference_v3', REFERENCE_PATH)
_reference = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_reference)


def _integer(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError(name+' must be an integer, excluding bool.')
    return int(value)


def _exact_rational(value, name):
    if type(value) is Q:
        # Fractions can retain NumPy scalar components. Canonicalize before
        # subtraction, exponentiation, or reference-transform arithmetic.
        return Q(int(value.numerator), int(value.denominator))
    return Q(_integer(value, name))


def _tensor(values, name):
    a = np.asarray(values, dtype=object)
    if not 1 <= a.ndim <= 3 or any(n < 1 or n > 257 for n in a.shape):
        raise ValueError(name+' must have one to three nonempty axes of degree at most 256.')
    if a.size > 2**22:
        raise ValueError(name+' exceeds the explicit resource guard.')
    # Convert NumPy scalar integers BEFORE arithmetic: object dtype alone
    # does not prevent overflow within a retained np.int64 scalar.
    return np.fromiter((_integer(v, name) for v in a.flat), dtype=object,
                       count=a.size).reshape(a.shape)


def _ceil_div(n, d):
    return -((-n)//d)


def _round_ratio(values, numerator_scale, denominator, upper):
    fn = _ceil_div if upper else lambda n, d: n//d
    return np.fromiter((fn(int(v)*numerator_scale, denominator) for v in values.flat),
                       dtype=object, count=values.size).reshape(values.shape)


@lru_cache(maxsize=32)
def _rounded_transform(n, left, width, precision):
    T, denominator = _reference.transform(n, left, width)
    if type(denominator) is not int or denominator <= 0:
        raise ValueError('Invalid exact transform denominator.')
    if T.shape != (n+1, n+1) or any(int(v) < 0 for v in T.flat):
        raise ValueError('The exact axis transform must be nonnegative.')
    S = 1 << precision
    low = _round_ratio(T, S, denominator, False)
    high = _round_ratio(T, S, denominator, True)
    if any(not (int(lo)*denominator <= int(t)*S <= int(hi)*denominator)
           for lo, t, hi in zip(low.flat, T.flat, high.flat)):
        raise ArithmeticError('Transform enclosure failed.')
    low.setflags(write=False)
    high.setflags(write=False)
    return low, high


def _positive_interval_product(TL, TU, CL, CU, scale, column_block):
    """Enclose T*C at the SAME fixed binary scale.

    TL/scale <= T <= TU/scale, 0 <= TL, and CL/scale <= C <= CU/scale.
    Since T is nonnegative, lower and upper endpoints may be propagated
    separately. The sign of each endpoint selects the adverse T bound.
    """
    lower = _reference.integer_matmul(TL, np.maximum(CL, 0), column_block)
    lower += _reference.integer_matmul(TU, np.minimum(CL, 0), column_block)
    lower = _round_ratio(lower, 1, scale, False)
    upper = _reference.integer_matmul(TU, np.maximum(CU, 0), column_block)
    upper += _reference.integer_matmul(TL, np.minimum(CU, 0), column_block)
    upper = _round_ratio(upper, 1, scale, True)
    if any(int(lo) > int(hi) for lo, hi in zip(lower.flat, upper.flat)):
        raise ArithmeticError('Invalid propagated interval.')
    return lower, upper


def to_bernstein_interval(C, denominator, bounds, *, precision=384,
                          column_block=1024, progress=None):
    """Return (lower, upper, scale), enclosing every exact coefficient.

    The returned arrays have object dtype and Python-int entries. The scale
    is 2**precision throughout all axis products. Input arrays are copied.
    A low precision may return wide intervals; it never changes validity.
    """
    precision = _integer(precision, 'precision')
    if not 1 <= precision <= 4096:
        raise ValueError('precision must lie in [1,4096].')
    column_block = _integer(column_block, 'column_block')
    if not 1 <= column_block <= 65536:
        raise ValueError('column_block must lie in [1,65536].')
    denominator = _integer(denominator, 'denominator')
    if denominator <= 0:
        raise ValueError('denominator must be positive.')
    C = _tensor(C, 'C')
    if not isinstance(bounds, (list, tuple)) or len(bounds) != C.ndim:
        raise ValueError('One exact interval per tensor axis is required.')
    parsed = []
    for b in bounds:
        if not isinstance(b, (list, tuple)) or len(b) != 2:
            raise ValueError('Each bound must be a pair.')
        left, right = (_exact_rational(v, 'bound') for v in b)
        if not 0 <= left < right:
            raise ValueError('Every axis requires 0 <= left < right.')
        parsed.append((left, right))
    if progress is not None and not callable(progress):
        raise TypeError('progress must be callable or None.')
    scale = 1 << precision
    low = _round_ratio(C, scale, denominator, False)
    high = _round_ratio(C, scale, denominator, True)
    del C
    for axis, (left, right) in enumerate(parsed):
        started = time.monotonic()
        n = low.shape[axis]-1
        TL, TU = _rounded_transform(n, left, right-left, precision)
        moved_low, moved_high = np.moveaxis(low, axis, 0), np.moveaxis(high, axis, 0)
        shape = moved_low.shape
        new_low, new_high = _positive_interval_product(
            TL, TU, moved_low.reshape(n+1, -1), moved_high.reshape(n+1, -1),
            scale, column_block)
        low = np.moveaxis(new_low.reshape(shape), 0, axis)
        high = np.moveaxis(new_high.reshape(shape), 0, axis)
        if progress is not None:
            progress(dict(axis=axis, degree=n, seconds=time.monotonic()-started,
                          precision=precision,
                          maximum_endpoint_bits=max(abs(int(v)).bit_length()
                                                    for a in (low, high) for v in a.flat),
                          maximum_width_units=max(int(hi)-int(lo)
                                                  for lo, hi in zip(low.flat, high.flat))))
    return low, high, scale


def to_int64_enclosure(low, high, scale, *, maximum_precision=50):
    """Outward-round an exact interval tensor to the old subdivision format.

    Returns (int64_lower, int64_upper, 2**p, p). Every integer magnitude is
    checked to be strictly below 2**59 BEFORE the int64 conversion. If no
    nonnegative output precision fits, raises OverflowError, never clips.
    """
    low, high = _tensor(low, 'lower'), _tensor(high, 'upper')
    if low.shape != high.shape or any(int(a) > int(b) for a, b in zip(low.flat, high.flat)):
        raise ValueError('Mismatched or inverted interval endpoints.')
    scale = _integer(scale, 'scale')
    if scale <= 0 or scale & (scale-1):
        raise ValueError('scale must be a positive power of two.')
    maximum_precision = _integer(maximum_precision, 'maximum_precision')
    if not 0 <= maximum_precision <= 50:
        raise ValueError('maximum_precision must lie in [0,50].')
    lo_min, hi_max = min(map(int, low.flat)), max(map(int, high.flat))
    limit = 1 << 59
    for precision in range(maximum_precision, -1, -1):
        output_scale = 1 << precision
        L, U = lo_min*output_scale//scale, _ceil_div(hi_max*output_scale, scale)
        if max(abs(L), abs(U)) < limit:
            lo = _round_ratio(low, output_scale, scale, False)
            hi = _round_ratio(high, output_scale, scale, True)
            if not all(-limit < int(v) < limit for a in (lo, hi) for v in a.flat):
                raise OverflowError('Final integer overflow guard failed.')
            return lo.astype(np.int64), hi.astype(np.int64), output_scale, precision
    raise OverflowError('Coefficient dynamic range exceeds the int64 handoff.')
