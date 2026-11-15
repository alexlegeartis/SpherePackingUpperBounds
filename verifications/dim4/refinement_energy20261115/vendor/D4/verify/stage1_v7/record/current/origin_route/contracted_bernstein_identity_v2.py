"""Exact Bernstein evaluation by contracting the transform before the tensor.

This validates the basis-transform formula at prescribed rational weights.
It does not inspect, validate or call an interval Bernstein converter.
No sign decision or global-certificate status is returned by this module.
"""
from fractions import Fraction as Q
from functools import lru_cache
from numbers import Integral, Rational
import math
import numpy as np


def _rational(value):
    assert isinstance(value, Rational) and not isinstance(value, (bool, float, np.floating))
    value = Q(value)
    # Fraction may retain native NumPy integer components. Normalize both
    # components before powers, products, subtraction or cached arithmetic.
    return Q(int(value.numerator), int(value.denominator))


def _integer(value):
    assert isinstance(value, Integral) and not isinstance(value, (bool, np.bool_))
    return int(value)


@lru_cache(maxsize=256)
def _axis_weights(degree, left, right, parameter):
    """Numerators of T^T b, over one positive common denominator.

    T[i,k] is the degree-n Bernstein coefficient of (left+width*t)^k.
    b[i] is B_i^n(parameter).  First contract b against binom(i,j),
    then apply the monomial translation/elevation coefficients.  All sums
    are exact integers; only a final common gcd reduction is performed.
    """
    n = degree
    p, q = parameter.numerator, parameter.denominator
    basis = tuple(math.comb(n, i) * p**i * (q-p)**(n-i)
                  for i in range(n+1))
    moments = tuple(sum(basis[i] * math.comb(i, j)
                        for i in range(j, n+1)) for j in range(n+1))
    width = right-left
    grid = math.lcm(left.denominator, width.denominator)
    a, w = int(left*grid), int(width*grid)
    common = math.lcm(*(math.comb(n, j) for j in range(n+1)))
    numerators = tuple(grid**(n-k) * sum(
        math.comb(k, j) * a**(k-j) * w**j
        * (common//math.comb(n, j)) * moments[j]
        for j in range(k+1)) for k in range(n+1))
    denominator = grid**n * common * q**n
    divisor = math.gcd(denominator, *numerators)
    assert denominator > 0 and divisor > 0
    return tuple(v//divisor for v in numerators), denominator//divisor


def contracted_axis_weights(degree, left, right, parameter):
    """Return exact T^T b as an integer vector and positive denominator.

    Validate before the cache lookup, so a float/bool equal to a previous
    rational cache key cannot bypass the exact-input contract.
    """
    degree = _integer(degree)
    left, right, parameter = map(_rational, (left, right, parameter))
    assert degree >= 0 and left < right and 0 <= parameter <= 1
    return _axis_weights(degree, left, right, parameter)


def clear_cache():
    _axis_weights.cache_clear()


def contracted_bernstein_evaluation(coefficients, denominator, bounds, parameters):
    """Evaluate (tensor T) C against tensor Bernstein weights exactly.

    The input C is an integer power-basis tensor divided by denominator.
    The returned Fraction is the contraction of its exact Bernstein tensor,
    computed without constructing that tensor.  Inputs are never mutated.
    """
    assert isinstance(coefficients, np.ndarray)
    assert coefficients.ndim >= 1 and all(n > 0 for n in coefficients.shape)
    denominator = _integer(denominator)
    assert denominator > 0
    assert len(bounds) == len(parameters) == coefficients.ndim
    assert all(len(interval) == 2 for interval in bounds)
    values = []
    python_integers = coefficients.dtype == object
    for value in coefficients.flat:
        number = _integer(value)
        python_integers = python_integers and type(value) is int
        values.append(number)
    # Native NumPy int64 arithmetic must never enter an exact contraction.
    # The usual object/Python-int input can be used directly; other integer
    # input dtypes are copied into canonical Python integers first.
    tensor = coefficients if python_integers else np.array(values, dtype=object).reshape(coefficients.shape)
    del values
    weights = []
    for n, interval, parameter in zip(coefficients.shape, bounds, parameters):
        vector, divisor = contracted_axis_weights(n-1, interval[0], interval[1], parameter)
        weights.append(np.array(vector, dtype=object))
        denominator *= divisor
    answer = tensor
    for vector in reversed(weights):
        answer = np.tensordot(answer, vector, axes=(-1, 0))
    return Q(int(answer), denominator)
