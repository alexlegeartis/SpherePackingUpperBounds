import numpy as np,json
D=np.load('/tmp/ct36/data.npz');H=D['H'];V=np.load('/tmp/ct36/min6.npy');A=V.T@V@H
assert np.array_equal(A@A@A,288*A@A-18432*A)
# Common denominator18432: projector0=I-3A/192 +A²/18432; P1=2A/96-A²/9216; P2=A²/18432-A/192
parts=np.array([H@(18432*np.eye(36,dtype=np.int64)-288*A+A@A),H@(384*A-2*A@A),H@(A@A-96*A)])
assert np.array_equal(parts.sum(0),18432*H)
np.save('/tmp/ct36/parts.npy',parts)
q=np.stack([np.einsum('ij,jk,ik->i',V,F,V)for F in parts],1);print('min6 types',np.unique(q,axis=0)/18432)
