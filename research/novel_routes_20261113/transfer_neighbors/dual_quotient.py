"""Find short explicit witnesses for every class of L*/L.

No complete enumeration is required when every class receives a witness.
The minimum6 of L is a catalogue input. Covering of L*/L by norm<=4 vectors
is a finite exact claim, independently replayable from the saved array.
"""
from pathlib import Path
import argparse,json,subprocess,time
import numpy as np,sympy as sp
from fpylll import IntegerMatrix,GSO,LLL,BKZ
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
ap=argparse.ArgumentParser();ap.add_argument('index',type=int);ap.add_argument('--pairs',type=int,default=2000000);args=ap.parse_args();i=args.index
G=np.loadtxt(P/f'Q32_{i}.txt',dtype=np.int64);D=2*sp.Matrix(G).inv()
assert all(x.q==1 for x in D);D=np.array(D,dtype=np.int64)
Q=IntegerMatrix.from_matrix(D.tolist());U=IntegerMatrix.identity(32);M=GSO.Mat(Q,gram=True,U=U,update=True);lll=LLL.Reduction(M);lll();BKZ.Reduction(M,lll,BKZ.Param(block_size=22,max_loops=2))()
H=np.array([[Q[j,k]for k in range(32)]for j in range(32)],dtype=np.int64);T=np.array([[U[j,k]for k in range(32)]for j in range(32)],dtype=np.int64)
assert np.array_equal(T@D@T.T,H);assert abs(sp.Matrix(T).det())==1
C=T@D # 2 times the dual-reduced basis expressed in original primal basis.
A=C.copy()%2;r=0;piv=[]
for j in range(32):
 v=np.flatnonzero(A[r:,j])
 if not len(v):continue
 t=r+v[0];A[[r,t]]=A[[t,r]]
 for k in range(32):
  if k!=r and A[k,j]:A[k]^=A[r]
 piv.append(j);r+=1
 if r==32:break
assert r==16
tags=np.array([sum(int(C[j,k]%2)<<a for a,k in enumerate(piv))for j in range(32)],dtype=np.uint16)
np.savetxt(P/f'dual_{i}_Q.txt',H,fmt='%d');np.savetxt(P/f'dual_{i}_R.txt',np.linalg.cholesky(H).T,fmt='%.18g')
pool=P/f'dual_minvecs_{i}.txt'
if not pool.exists():
 q=subprocess.run([str(ROOT/'research/density_breakthroughs_20261113/bin/shell_discovery'),'32',str(P/f'dual_{i}_Q.txt'),str(P/f'dual_{i}_R.txt'),'6','200000000',str(pool),'0'],capture_output=True,text=True,check=True)
 (P/f'dual_enum_{i}.json').write_text(q.stdout)
raw=np.loadtxt(pool,dtype=np.int64);V=raw[:,1:];norm=raw[:,0];assert np.all(np.einsum('bi,ij,bj->b',V,H,V)==norm)
assert np.max(norm)<=6
labels=np.zeros(len(V),dtype=np.uint16)
for j in range(32):labels^=np.where((V[:,j]&1)!=0,tags[j],np.uint16(0))
W=np.zeros((65536,32),dtype=np.int16);covered=np.zeros(65536,dtype=bool);covered[0]=True
labs,ix=np.unique(labels,return_index=True);W[labs]=V[ix];covered[labs]=True
start=time.time();rng=np.random.default_rng(7200+i);VH=V@H;tested=0
for step in range((args.pairs+99999)//100000):
 a=rng.integers(0,len(V),100000);b=rng.integers(0,len(V),100000);ip=np.einsum('bi,bi->b',V[a],VH[b]);q=norm[a]+norm[b]-2*np.abs(ip);sg=np.where(ip>=0,1,-1)
 ll=labels[a]^labels[b];mask=(q<=8)&(~covered[ll]);llnew,ixnew=np.unique(ll[mask],return_index=True);ids=np.flatnonzero(mask)[ixnew]
 W[llnew]=(V[a[ids]]-sg[ids,None]*V[b[ids]]).astype(np.int16);covered[llnew]=True;tested+=100000
 if step%5==0:print('pairs',tested,'covered',int(covered.sum()),flush=True)
 if covered.all():break
# Directed pairing: seek a short sum in one currently uncovered class.
bylabel={int(k):[] for k in labs}
for j,k in enumerate(labels):bylabel[int(k)].append(j)
target_tests=0
for target in np.flatnonzero(~covered):
 for attempt in range(400):
  a=int(rng.integers(len(V)));want=int(target)^int(labels[a]);choices=bylabel.get(want)
  if not choices:continue
  b=choices[int(rng.integers(len(choices)))];ip=int(V[a]@VH[b]);q=int(norm[a]+norm[b]-2*abs(ip));target_tests+=1
  if q<=8:
   W[target]=V[a]-(1 if ip>=0 else -1)*V[b];covered[target]=True;break
missing=np.flatnonzero(~covered).tolist();out=P/f'dual_quotient_{i}.npz'
assert np.max(np.abs(W))<128
np.savez_compressed(out,G=G,D=D,H=H,T=T,C=C,piv=np.array(piv),witnesses=W.astype(np.int8),covered=covered)
report={'index':i,'classes':65536,'covered':int(covered.sum()),'missing':missing,'pool_size':len(V),'random_pairs':tested,'directed_pairs':target_tests,'maximum_witness_abs_coordinate':int(np.abs(W).max()),'seconds_after_enumeration':time.time()-start,'artifact':out.name,'scope':'Every covered class has an explicit exact norm<=4 dual witness; enumeration completeness is unnecessary for full coverage.'}
(P/f'dual_quotient_{i}.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
