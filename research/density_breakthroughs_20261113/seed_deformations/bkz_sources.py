import sys,pathlib,numpy as np,time,json
from fpylll import IntegerMatrix,BKZ
HERE=pathlib.Path(__file__).resolve().parent
for name in sys.argv[1:]:
 p=HERE/name;d=np.load(p/'seed.npz');G=d['G'];n=len(G)
 B0=np.rint(np.linalg.cholesky(G)*1e8).astype(np.int64)
 B=IntegerMatrix.from_matrix(B0.tolist());U=IntegerMatrix.identity(n)
 t=time.time();BKZ.reduction(B,BKZ.Param(block_size=22,max_loops=2),U=U,float_type='long double')
 U=np.array([[int(U[i,j]) for j in range(n)]for i in range(n)],np.int64);H=U@G@U.T
 assert np.array_equal(U@B0,np.array([[int(B[i,j])for j in range(n)]for i in range(n)]))
 print(name,'time',time.time()-t,'diagonal',H.diagonal(),flush=True)
 np.savez(p/'reduced.npz',G=G,H=H,U=U)
 np.savetxt(p/'Qred.txt',H,fmt='%d');np.savetxt(p/'Rred.txt',np.linalg.cholesky(H).T,fmt='%.17g')
