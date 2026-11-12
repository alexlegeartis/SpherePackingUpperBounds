"""Rebuild large section contact arrays from the shared discovery shell."""
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
raw=np.loadtxt(P.parent/'contact_walk/contacts_reduced.txt',dtype=np.int64)
U=np.loadtxt(ROOT/'constructions/dim36/data/U.txt',dtype=np.int64)
Q=np.loadtxt(ROOT/'constructions/dim36/data/Q_original.txt',dtype=np.int64)
G4=np.loadtxt(ROOT/'constructions/dim36/data/G4.txt',dtype=np.int64)
V=raw[:,1:]@U
assert np.all(raw[:,0]==160) and np.all(np.einsum('ij,ij->i',V@Q,V)==160)
oldnorm=np.einsum('ij,ij->i',V@G4,V)
np.savez_compressed(P/'full_contacts.npz',V=V,oldnorm=oldnorm)
D=np.load(P/'dual_lll.npz')
np.savez_compressed(P/'normal_contact_coordinates.npz',V=V@D['U'].T,H=D['H'])
print('exact individual norm-eight contact representatives',len(V))
