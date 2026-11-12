import numpy as np,pathlib,json,sympy as sp
from scipy.linalg import qr,svd
from fractions import Fraction
p=pathlib.Path(__file__).resolve().parent/'dim48mod3';H=np.load(p/'reduced.npz')['H'];V=np.load(p/'contacts.npy');n=len(H);rng=np.random.default_rng(42);rng.shuffle(V);i,j=np.triu_indices(n);M=(V[:6000,i]*V[:6000,j]*(1+(i!=j))).astype(float)
u,s,vh=svd(M,full_matrices=False,lapack_driver='gesdd');N=vh[s<1e-7].T;print('nullity',N.shape[1],flush=True)
_,_,piv=qr(N.T,pivoting=True,mode='economic');free=piv[:N.shape[1]];N=N@np.linalg.inv(N[free]);print('free',[(int(i[x]),int(j[x]))for x in free],flush=True)
Ds=[]
for k in range(N.shape[1]):
 vals=[Fraction(float(x)).limit_denominator(1024)for x in N[:,k]];den=int(sp.ilcm(*[x.denominator for x in vals]));vec=np.array([int(x*den)for x in vals],dtype=np.int64);D=np.zeros((n,n),np.int64);D[i,j]=vec;D[j,i]=vec
 res=np.einsum('ij,jk,ik->i',V,D,V);print(k,'den',den,'maxD',abs(D).max(),'resid',abs(res).max(),flush=True);assert not np.any(res);Ds.append(D)
Ds=np.array(Ds);np.savez(p/'flexes.npz',H=H,Ds=Ds,V=V,free=free)
print('traces',[np.trace(np.linalg.solve(H,D))for D in Ds],flush=True)
C=np.linalg.cholesky(H);Ci=np.linalg.inv(C);A=np.array([Ci@D@Ci.T for D in Ds]);print('eigenranges',[(np.linalg.eigvalsh(z)[0],np.linalg.eigvalsh(z)[-1])for z in A])
(p/'flex_report.json').write_text(json.dumps({'nullity':len(Ds),'contacts_checked':len(V),'integer_directions':Ds.tolist(),'free_entries':[(int(i[x]),int(j[x]))for x in free],'trace_directions':[float(np.trace(np.linalg.solve(H,D)))for D in Ds]},indent=2))
