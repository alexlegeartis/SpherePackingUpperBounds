#!/usr/bin/env python3
"""Discover and exactly certify coordinate-normal sections inside Gamma72's Leech block."""
from pathlib import Path
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
import json
P=Path(__file__).resolve().parent
G=np.load(P.parent/'higher_parent'/'G72.npy')
F=G[:24,:24]//2
assert np.array_equal(2*F,G[:24,:24])
assert sp.Matrix(F).det()==1
raw=np.loadtxt(P/'leech_minvecs.txt',dtype=np.int64)
V=raw[:,1:]
assert np.all(np.einsum('ij,jk,ik->i',V,F,V)==4)
VF=V@F
rng=np.random.default_rng(20261113)
target={1:4,2:12,3:32,4:64,5:128,6:192,7:256,8:256,9:512,10:768,11:1024,12:1024,13:1024,14:768,15:512,16:256}
best={}
for attempt in range(60):
 B=np.empty((0,24),dtype=np.int64)
 residual=np.full(len(V),4.0)
 Q=np.empty((0,24))
 for k in range(1,17):
  pos=residual>1e-7
  minimum=residual[pos].min()
  inds=np.flatnonzero(pos & (residual<minimum+1e-7))
  row=V[rng.choice(inds)]
  ortho=row.astype(float)-(Q@F@row)@Q if len(Q) else row.astype(float)
  ortho/=np.sqrt(ortho@F@ortho)
  Q=np.vstack((Q,ortho));residual-=(VF@ortho)**2
  B=np.vstack((B,row))
  K=B@F@B.T; det=int(sp.Matrix(K).det())
  if k not in best or det<best[k][0]:
   best[k]=(det,B.copy())
   print('improved',attempt,k,det,'target',target[k],flush=True)
  if det>target[k]: break
 if all(k in best and best[k][0]<=target[k] for k in target): break
out={}
for k,(d,B24) in sorted(best.items()):
 B=np.zeros((k,72),dtype=np.int64);B[:,:24]=B24
 K=B@G@B.T
 smith=smith_normal_form(sp.Matrix(B),domain=ZZ)
 factors=[abs(int(smith[i,i])) for i in range(k)]
 assert factors==[1]*k,(k,factors)
 assert int(sp.Matrix(K).det())==2**k*d
 n=72-k
 density_sq=sp.Rational(2)**n/sp.det(sp.Matrix(K))
 out[n]={'normal_rank':k,'leech_normal_det':d,'normal_det':int(sp.Matrix(K).det()),'smith_factors':factors,'center_density_sq':str(density_sq),'center_density':float(sp.sqrt(density_sq)),'inherited_minimum':8}
 np.savez_compressed(P/f'gamma72_leech_section_rank{k}.npz',G=G,B=B,K=K)
(P/'gamma72_leech_sections.json').write_text(json.dumps(out,indent=2)+'\n')
