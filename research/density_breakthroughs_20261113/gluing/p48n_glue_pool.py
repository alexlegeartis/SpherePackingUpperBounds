"""Exact P48n minimum-vector pool from saved E8 gluing data.
All coordinates are rows in the original P48n Gram basis; no completeness claim.
"""
from pathlib import Path
import numpy as np,sympy as sp,json,time
P=Path(__file__).resolve().parent;K=P.parents[3]/'KissingNumbers';S=K/'verifications/improved/dim41-46-p48-heads/population_replay/data';j=np.load(S/'prep_e8.npz');d=np.load(S/'e8_complement_reduced.npz');b=d['b'];q=d['q'];BB=j['BB'];G=j['G'];Ginv=np.array(sp.Matrix(G).inv(),dtype=np.int64);qi3=np.array(3*sp.Matrix(q).inv(),dtype=np.int64);assert np.array_equal(G@Ginv,np.eye(48,dtype=np.int64))
start=time.time();uv=[]
def go(pre,r,par):
 if len(pre)==8:
  if sum(pre)%4==0:uv.append(pre)
  return
 for x in range(-6+par,7,2):
  if x*x<=r:go(pre+[x],r-x*x,par)
go([],32,0);go([],32,1);u=np.array(uv,dtype=np.int64);un=np.sum(u*u,axis=1);assert len(u)==26641
pair=u@j['A'].T;assert not np.any(pair%4);lifts=(pair//4)@j['U'][:8]@BB;dots=lifts@b.T;assert not np.any(dots%96);c=(dots//96)@qi3;normal3=3*lifts-c@b
assert np.array_equal(np.sum(normal3*normal3,axis=1),72*un)
classes={}
for i in np.flatnonzero((un==24)|(un==32)):
 classes.setdefault((int(un[i]),(c[i]%3).astype(np.uint8).tobytes()),[]).append(i)
H=np.load(S/'dual40_min.npz')['H'].astype(np.int64);hprod=H@b.T;assert not np.any(hprod%96);C=(hprod//96)@qi3;assert not np.any(C%3);C//=3;assert np.array_equal(C@b,H);assert np.all(np.einsum('bi,ij,bj->b',C,q,C)==30)
AUX=np.load(K/'verifications/improved/dim40-affine-shells/source/auxiliary.npz')['norm36'].astype(np.int64);qold=np.array(json.loads((K/'verifications/improved/dim40-affine-shells/source/gram.json').read_text()));assert np.array_equal(qold,q)
parts=[j['B'],np.load(P/'contacts40.npz')['V']@b];counts={}
for name,coeff,norm in [('short',C,32),('middle',AUX,24)]:
 out=[]
 for v in coeff:
  inds=classes.get((norm,(v%3).astype(np.uint8).tobytes()),[])
  for i in inds:
   p3=v@b+normal3[i];assert not np.any(p3%3);out.append(p3//3)
 out=np.array(out,dtype=np.int64);assert np.all(np.sum(out*out,axis=1)==576);parts.append(out);counts[name]=len(out);print(name,len(out),'seconds',time.time()-start,flush=True)
X=np.unique(np.concatenate(parts),axis=0);assert np.all(np.sum(X*X,axis=1)==576);znum=(X@BB.T)@Ginv;assert not np.any(znum%96);V=znum//96;sign=np.sign(V[np.arange(len(V)),np.argmax(V!=0,axis=1)]);V=np.unique(V*sign[:,None],axis=0);assert np.all(np.einsum('bi,ij,bj->b',V,G,V)==6)
np.savez_compressed(P/'P48n_glue_pool.npz',G=G,V=V,E8=j['T']);(P/'P48n_glue_pool.json').write_text(json.dumps({'points_up_to_sign':len(V),'source_groups':counts,'exact_norm':6,'minimum_complete':False,'seconds':time.time()-start},indent=2));print('DONE',len(V),'time',time.time()-start,flush=True)
