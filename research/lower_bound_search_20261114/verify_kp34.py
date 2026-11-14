"""Exact minimum, perfection and positive eutaxy of a cyclic-code 34D seed.

The contact set is generated algebraically, without a shortest-vector search.
This certifies a strict local density maximum among lattices. It does not
certify a sphere-packing record or exclude other packings.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import numpy as np
from verify_codes import MUL

P = Path(__file__).resolve().parent


def determinant(matrix):
    A = [[int(x) for x in row] for row in matrix]
    sign, previous = 1, 1
    for k in range(len(A)-1):
        if not A[k][k]:
            pivot = next(i for i in range(k+1, len(A)) if A[i][k])
            A[k], A[pivot] = A[pivot], A[k]
            sign = -sign
        diagonal = A[k][k]
        for i in range(k+1, len(A)):
            for j in range(k+1, len(A)):
                numerator = diagonal*A[i][j]-A[i][k]*A[k][j]
                assert numerator % previous == 0
                A[i][j] = numerator // previous
            A[i][k] = 0
        previous = diagonal
    return sign*A[-1][-1]


def code():
    g = np.array([1, 2, 0, 2, 2, 2, 0, 2, 1], dtype=np.uint8)
    even = np.r_[g, 0] ^ np.r_[0, g]
    G = np.array([np.pad(even, (i, 7-i)) for i in range(8)], dtype=np.uint8)
    C = np.zeros((1, 17), dtype=np.uint8)
    for row in G:
        C = np.concatenate([C ^ MUL[a, row] for a in range(4)])
    weight = np.count_nonzero(C, axis=1)
    assert len(C) == 4**8 and np.count_nonzero(weight == 0) == 1
    assert weight[1:].min() == 8 and np.count_nonzero(weight == 8) == 1530
    return C, weight


def generate_contacts(B, C, weight):
    # Each nonzero quaternary symbol and nonzero trit has one norm-one lift.
    unit = np.zeros((4, 3, 2), dtype=np.int64)
    unit[1, 1], unit[1, 2] = (1, 0), (-1, 0)
    unit[2, 1], unit[2, 2] = (0, 1), (0, -1)
    unit[3, 1], unit[3, 2] = (-1, -1), (1, 1)
    patterns = np.array(list(itertools.product((1, 2), repeat=8)), dtype=np.int64)
    patterns = patterns[patterns.sum(axis=1) % 3 == 0]
    assert len(patterns) == 86
    blocks = []
    for word in C[weight == 8]:
        positions = np.flatnonzero(word)
        W = np.zeros((86, 17, 2), dtype=np.int64)
        for index, position in enumerate(positions):
            W[:, position] = unit[int(word[position]), patterns[:, index]]
        blocks.append(W.reshape(86, 34))
    zero_q = []
    for left, right in itertools.combinations(range(17), 2):
        for trit in (1, 2):
            for a, b in itertools.product((1, 2, 3), repeat=2):
                W = np.zeros((17, 2), dtype=np.int64)
                # Doubling a unit with trit -t gives a zero-q norm-four lift
                # with trit t, because 2*(-t) == t modulo three.
                W[left] = 2*unit[a, 3-trit]
                W[right] = 2*unit[b, trit]
                zero_q.append(W.reshape(34))
    assert len(zero_q) == 2448
    Y = np.concatenate(blocks + [np.array(zero_q)])
    assert len(Y) == 134028 and not (Y.sum(axis=1) % 3).any()
    # The supplied inverse is an integer witness, checked before use.
    inverse6 = np.loadtxt(P / 'inputs/kp34_inverse6.txt', dtype=np.int64)
    assert np.array_equal(B @ inverse6, 6*np.eye(34, dtype=np.int64))
    numerator = Y @ inverse6
    assert not (numerator % 6).any()
    V = numerator // 6
    assert np.array_equal(V @ B, Y)
    return V


def rank_mod_prime(A, prime):
    A = A.copy() % prime
    rank = 0
    for column in range(A.shape[1]):
        choices = np.flatnonzero(A[rank:, column])
        if not len(choices):
            continue
        pivot = rank + int(choices[0])
        A[[rank, pivot]] = A[[pivot, rank]]
        multipliers = A[rank+1:, column] * pow(int(A[rank, column]), -1, prime) % prime
        A[rank+1:] = (A[rank+1:] - multipliers[:, None]*A[rank]) % prime
        rank += 1
        if rank == A.shape[0]:
            break
    return rank


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=P / 'work/kp34.json')
    args = parser.parse_args()
    B = np.loadtxt(P / 'inputs/kp34_basis.txt', dtype=np.int64)
    assert B.shape == (34, 34)
    C, weight = code()
    encode_weights = 1 << (2*np.arange(17, dtype=np.uint64))
    encoded = set(map(int, C.astype(np.uint64) @ encode_weights))
    q = (B[:, 0::2] % 2) + 2*(B[:, 1::2] % 2)
    assert all(int(value) in encoded for value in q.astype(np.uint64) @ encode_weights)
    assert not (B.sum(axis=1) % 3).any()
    det_B = abs(determinant(B))
    assert det_B == 3*4**9
    metric2 = np.kron(np.eye(17, dtype=np.int64), np.array([[2, -1], [-1, 2]], dtype=np.int64))
    twice_Q = B @ metric2 @ B.T
    assert not (twice_Q % 2).any()
    Q = twice_Q // 2
    det_Q = determinant(Q)
    assert det_Q == 4649045868
    V = generate_contacts(B, C, weight)
    assert np.all(np.einsum('vi,ij,vj->v', V, Q, V) == 8)
    # A nonzero q has at least eight nonzero coordinates, hence norm >= 8.
    # For q=0, a vector of norm below eight could only be a doubled unit
    # at one coordinate. Its trit is nonzero, contradicting parity zero.
    # Thus the algebraic code construction itself proves minimum eight.
    moment = V.T @ V @ Q
    assert np.array_equal(moment, 31536*np.eye(34, dtype=np.int64))
    selected = np.loadtxt(P / 'inputs/kp34_perfection_rows.txt', dtype=np.int64)
    assert selected.shape == (595,) and len(set(map(int, selected))) == 595
    W = V[selected]
    i, j = np.triu_indices(34)
    A = W[:, i]*W[:, j]*np.where(i == j, 1, 2)
    assert rank_mod_prime(A, 1000003) == 595
    result = {'dimension': 34, 'minimum_squared': 8, 'Gram_determinant': det_Q,
              'center_density_exact': '4^8/3^(19/2)', 'explicit_contacts': len(V),
              'positive_eutaxy_moment_scalar': 31536,
              'perfection_rank_mod_1000003': 595,
              'contact_sha256_le_i64': hashlib.sha256(V.astype('<i8').tobytes()).hexdigest(),
              'conclusion': 'Strict local density maximum among lattices; not a new record.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    print('PASS: exact minimum, positive eutaxy and perfection of the cyclic-code 34D seed.')


if __name__ == '__main__':
    main()
