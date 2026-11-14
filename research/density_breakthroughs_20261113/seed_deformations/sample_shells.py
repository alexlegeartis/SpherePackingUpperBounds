import numpy as np,pathlib,subprocess,json,sys
from fpylll import IntegerMatrix,BKZ
HERE=pathlib.Path(__file__).resolve().parent
name=sys.argv[1];trials=int(sys.argv[2]);ncap=int(sys.argv[3]);bound=float(sys.argv[4]);p=HERE/name
G=np.load(p/'reduced.npz')['H'];n=len(G);gen=np.random.default_rng(13)
V=np.loadtxt(p/'shellred.txt',dtype=np.int64,ndmin=2)[:,1:]
allv=set()
def add(Z):
 for v in Z:
  k=np.flatnonzero(v)[0]
  if v[k]<0:v=-v
  allv.add(tuple(map(int,v)))
add(V);recs=[]
for trial in range(trials):
 C=np.linalg.cholesky(G);perm=gen.permutation(n);B0=np.rint(C[perm]*1e8).astype(np.int64);B=IntegerMatrix.from_matrix(B0.tolist());T=IntegerMatrix.identity(n)
 BKZ.reduction(B,BKZ.Param(block_size=20,max_loops=2),U=T,float_type='long double')
 T=np.array([[int(T[i,j])for j in range(n)]for i in range(n)],dtype=np.int64)@np.eye(n,dtype=np.int64)[perm];Q=T@G@T.T
 np.savetxt(p/'Qsample.txt',Q,fmt='%d');np.savetxt(p/'Rsample.txt',np.linalg.cholesky(Q).T,fmt='%.17g')
 run=subprocess.run([str(pathlib.Path(__file__).resolve().parent.parent/'bin'/'shell_discovery'),str(n),str(p/'Qsample.txt'),str(p/'Rsample.txt'),str(bound),str(ncap),str(p/'shellsample.txt'),'0'],capture_output=True,text=True)
 Z=np.loadtxt(p/'shellsample.txt',dtype=np.int64,ndmin=2)[:,1:]@T
 assert np.all(np.einsum('ij,jk,ik->i',Z,G,Z)==int(bound))
 add(Z);r=json.loads(run.stdout);r['trial']=trial;r['total']=len(allv);recs.append(r);print(json.dumps(r),flush=True)
 np.save(p/'contacts.npy',np.array(sorted(allv),dtype=np.int64));(p/'sampling.json').write_text(json.dumps(recs,indent=2))
