#!/usr/bin/env python3
"""Exact eutaxy and perfection certificate for the existing 36D lattice.

This is not a new density record. It uses explicit norm-160 contacts of Q;
the existing construction's minimum certificate supplies min(Q)=160.
Neither the contact identity nor the modular rank test requires a complete
enumeration of the minimal shell.
"""
from pathlib import Path
import json
import numpy as np

P=Path(__file__).resolve().parent
S=P.parents[2]/'constructions/dim36'
Q=np.loadtxt(S/'Q_reduced.txt',dtype=np.int64)
U=np.loadtxt(S/'data/U.txt',dtype=np.int64)
D=U@np.loadtxt(S/'data/D.txt',dtype=np.int64)@U.T
Z=np.load(P/'contacts36.npz')['Z'].astype(np.int64)
assert np.all(np.einsum('bi,ij,bj->b',Z,Q,Z)==160)
slopes=np.einsum('bi,ij,bj->b',Z,D,Z)
den=417590880
coefficients={-8:1326,-4:2581,-2:698,0:1326}
assert set(map(int,slopes))==set(coefficients)
M=np.zeros((36,36),dtype=np.int64)
for slope,c in coefficients.items():
    V=Z[slopes==slope];M+=c*(V.T@V)
assert 36*int(abs(M).max())*int(abs(Q).max())<2**63
assert np.array_equal(M@Q,den*np.eye(36,dtype=np.int64))

prime=65521
assert all(prime%d for d in range(2,int(prime**.5)+1))
certificate=P/'perfection36.json'
if certificate.exists():
    row_indices=np.array(json.loads(certificate.read_text())['contact_indices'],dtype=int)
else:
    row_indices=np.random.default_rng(360013).choice(len(Z),2500,replace=False)
i,j=np.triu_indices(36)
B=(Z[row_indices,i[:,None]].T*Z[row_indices,j[:,None]].T)%prime
assert B.shape==(len(row_indices),666)
original=row_indices.copy();rank=0;chosen=[]
for c in range(666):
    candidates=np.flatnonzero(B[rank:,c])
    if not len(candidates):continue
    r=rank+int(candidates[0]);B[[rank,r]]=B[[r,rank]];original[[rank,r]]=original[[r,rank]]
    inv=pow(int(B[rank,c]),-1,prime)
    B[rank,c:]=(B[rank,c:]*inv)%prime
    B[rank+1:,c:]=(B[rank+1:,c:]-B[rank+1:,c,None]*B[rank,c:])%prime
    chosen.append(int(original[rank]));rank+=1
    if rank==len(B):break
assert rank==666
certificate.write_text(json.dumps({'prime':prime,'rank':rank,'contact_indices':chosen,'scope':'nonzero 666-by-666 contact determinant modulo prime proves full rank over the rationals'},indent=2))
report={'dimension':36,'exact_contacts':len(Z),'minimum_contact_norm_Q':160,'eutaxy_common_denominator':den,'positive_weight_numerators':coefficients,'exact_identity':'(sum_v numerator(v) vv^T) Q = 417590880 I','perfection_rank':rank,'rank_prime':prime,'conclusion':'perfect and eutactic using explicit contacts; with the prior minimum certificate, the existing 36D lattice is a strict local maximum among lattice packings modulo scale and isometry','new_density_record':False}
(P/'local36_verified.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
