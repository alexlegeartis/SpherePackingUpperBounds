"""Exact integer Bernstein transforms using FLINT matrix arithmetic.

For x=A/D+(B/D)t and a polynomial of degree n, the transform entry is
sum_j binom(k,j) binom(i,j)/binom(n,j) * (A/D)^(k-j) (B/D)^j.
With H=lcm_j binom(n,j), the denominator D^n H clears every term.
The formula below evaluates those integers directly, then cancels their
common gcd. It changes neither polynomial nor any inequality criterion.
FLINT multiplication is over Z; no floating point or modular guess enters.
"""
from fractions import Fraction as Q
from functools import lru_cache
import math
import numpy as np
from flint import fmpz_mat


@lru_cache(maxsize=64)
def transform(n, left, width):
    n = int(n)
    left, width = Q(left), Q(width)
    assert n >= 0 and width > 0
    D = math.lcm(left.denominator, width.denominator)
    A, B = int(left*D), int(width*D)
    choose = [[math.comb(i, j) for j in range(i+1)] for i in range(n+1)]
    H = math.lcm(*choose[n])
    pa, pb, pd = ([value**j for j in range(n+1)] for value in [A, B, D])
    h = [H//v for v in choose[n]]
    rows = [[pd[n-k]*sum(choose[k][j]*choose[i][j]*h[j]*pa[k-j]*pb[j]
                         for j in range(min(i, k)+1))
             for k in range(n+1)] for i in range(n+1)]
    denominator = pd[n]*H
    divisor = math.gcd(denominator, *(v for row in rows for v in row))
    assert divisor > 0
    matrix = np.array([[v//divisor for v in row] for row in rows], dtype=object)
    return matrix, denominator//divisor


def integer_matmul(left, right, column_block=1024):
    left, right = np.asarray(left, dtype=object), np.asarray(right, dtype=object)
    assert left.ndim == right.ndim == 2 and left.shape[1] == right.shape[0]
    assert isinstance(column_block, int) and column_block > 0
    assert all(isinstance(v, (int, np.integer)) for v in left.flat)
    assert all(isinstance(v, (int, np.integer)) for v in right.flat)
    n, k = left.shape
    m = right.shape[1]
    lmatrix = fmpz_mat(n, k, [int(v) for v in left.flat])
    result = np.empty((n, m), dtype=object)
    for start in range(0, m, column_block):
        end = min(m, start+column_block)
        rmatrix = fmpz_mat(k, end-start, [int(v) for v in right[:, start:end].flat])
        product = lmatrix*rmatrix
        result[:, start:end] = np.array([int(v) for v in product.entries()], dtype=object).reshape(n, end-start)
    return result


def to_bernstein(C, denominator, bounds):
    C = np.asarray(C, dtype=object)
    assert C.ndim == len(bounds) and isinstance(denominator, int) and denominator > 0
    for axis, (left, right) in enumerate(bounds):
        n = C.shape[axis]-1
        T, td = transform(n, Q(left), Q(right)-Q(left))
        moved = np.moveaxis(C, axis, 0)
        shape = moved.shape
        moved = integer_matmul(T, moved.reshape(n+1, -1)).reshape(shape)
        C = np.moveaxis(moved, 0, axis)
        denominator *= td
    return C, denominator
