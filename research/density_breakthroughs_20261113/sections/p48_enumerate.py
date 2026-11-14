"""Capped floating discovery enumeration of P48 minimal vectors."""
from pathlib import Path
import numpy as np,sys,json,subprocess,time
from fpylll import IntegerMatrix,BKZ
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
ambient=sys.argv[1]if len(sys.argv)>1 else 'P48p';p=P/f'{ambient}_enumeration';p.mkdir(exist_ok=True)
f=np.load(P/f'{ambient}_pool.npz');G=f['G'];n=len(G)
B=IntegerMatrix.from_matrix(np.rint(np.linalg.cholesky(G)*1e8).astype(np.int64).tolist());U0=IntegerMatrix.identity(n)
BKZ.reduction(B,BKZ.Param(block_size=28,max_loops=3),U=U0)
U=np.array([[int(U0[i,j])for j in range(n)]for i in range(n)],dtype=np.int64);Q=U@G@U.T
np.savetxt(p/'Q.txt',Q,fmt='%d');np.savetxt(p/'R.txt',np.linalg.cholesky(Q).T,fmt='%.17g');np.savez(p/'basis.npz',U=U,G=G,Q=Q)
print('basis',Q.diagonal(),flush=True)
r=subprocess.run([str(ROOT.parent/'shell_discovery'),str(n),str(p/'Q.txt'),str(p/'R.txt'),'6','300000000',str(p/'vectors.txt'),'0'],capture_output=True,text=True)
print(r.stdout,r.stderr,flush=True);(p/'enumeration.log').write_text(r.stdout+r.stderr)
raw=np.loadtxt(p/'vectors.txt',dtype=np.int64).reshape(-1,n+1);V=raw[:,1:]@U
assert np.all(np.einsum('bi,ij,bj->b',V,G,V)==6)
np.savez_compressed(p/'pool.npz',G=G,V=V)
print('POOL',len(V),flush=True)
