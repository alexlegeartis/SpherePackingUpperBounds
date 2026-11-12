"""Add exact norm8/10 normal vectors; minimum vectors alone need not span good normals."""
from pathlib import Path
import numpy as np,json,time
P=Path(__file__).resolve().parent;j=np.load(P/'P48n_glue_pool.npz');G=j['G'];V=j['V'].astype(np.int16);VG=V.astype(np.int64)@G;assert np.max(abs(VG))<32768;VG=VG.astype(np.int16);rng=np.random.default_rng(3840);out={8:[],10:[],12:[]};limits={8:150000,10:120000,12:30000};start=time.time()
while any(sum(len(a)for a in out[k])<limits[k] for k in out):
 a=rng.integers(len(V),size=200000);b=rng.integers(len(V),size=200000);ip=np.einsum('bi,bi->b',VG[a].astype(np.int64),V[b].astype(np.int64));
 for norm in out:
  want=(12-norm)//2;inds=np.flatnonzero(abs(ip)==want);used=sum(len(x)for x in out[norm]);inds=inds[:max(0,limits[norm]-used)]
  if len(inds):out[norm].append(V[a[inds]]-np.where(ip[inds]>=0,1,-1)[:,None]*V[b[inds]])
 print({k:sum(len(a)for a in out[k])for k in out},flush=True)
W=np.concatenate([V]+[np.concatenate(out[k])for k in out]).astype(np.int16);sg=np.sign(W[np.arange(len(W)),np.argmax(W!=0,axis=1)]);W*=sg[:,None];W=np.unique(W,axis=0);norms=np.einsum('bi,ij,bj->b',W.astype(np.int64),G,W.astype(np.int64));assert set(norms)<={6,8,10,12};np.savez_compressed(P/'P48n_extended_normal_pool.npz',G=G,V=W,E8=j['E8'],norms=norms);r={'vectors':len(W),'norm_histogram':{str(k):int(np.sum(norms==k))for k in np.unique(norms)},'seconds':time.time()-start,'exact_lattice_vectors':True,'complete':False};(P/'P48n_extended_normal_pool.json').write_text(json.dumps(r,indent=2));print('DONE',r,flush=True)
