import numpy as np,sympy as sp,json,subprocess
from pathlib import Path
P=Path(__file__).resolve().parent;rng=np.random.default_rng(362026);Q=np.load(P/'candidate_1.npy');out=[]
for it in range(16):
 U0=np.eye(36,dtype=np.int64)[rng.permutation(36)]
 for k in range(20):
  i,j=rng.choice(36,2,False);U0[i]+=rng.choice([-1,1])*U0[j]
 H=U0@Q@U0.T;B=np.rint(np.linalg.cholesky(H/20)*1e7).astype(np.int64);_,U1=sp.Matrix(B).lll_transform(delta=sp.Rational(99,100));U=np.array(U1,dtype=np.int64)@U0;H=U@Q@U.T
 np.savetxt(P/'Qtemp.txt',H,fmt='%d');np.savetxt(P/'Rtemp.txt',np.linalg.cholesky(H/20).T)
 ans=subprocess.run([str(P/'short'),str(P/'Rtemp.txt'),str(P/'Qtemp.txt'),'15000000'],capture_output=True,text=True).stdout;print(it,ans.strip(),flush=True);out.append({'trial':it,'U':U.tolist(),'result':ans})
 json.dump(out,open(P/'random_bases.json','w'),indent=2)
 if 'found 1' in ans:break
