"""Convert saved exact 36D contacts into the current reduced basis."""
from pathlib import Path
import numpy as np
import sympy as sp

P = Path(__file__).resolve().parent
SOURCE = P.parents[2] / 'constructions/dim36'
Q = np.loadtxt(SOURCE/'Q_reduced.txt', dtype=np.int64)
U = np.loadtxt(SOURCE/'data/U.txt', dtype=np.int64)
inverse = sp.Matrix(U).inv()
assert all(x.q == 1 for x in inverse)
inverse = np.array(inverse, dtype=np.int64)
V = np.load(SOURCE/'provenance/ct36_contact_walk_20261112/flex.npz')['V'] @ inverse
norms = np.einsum('bi,ij,bj->b', V, Q, V)
assert np.all(norms == 160)
np.savez_compressed(P/'contacts_reduced.npz', Q=Q, V=V, norms=norms)
print(f'Checked {len(V)} exact contacts in the reduced Gram basis.')
