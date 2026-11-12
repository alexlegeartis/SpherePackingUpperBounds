import numpy as np,pathlib,json,subprocess,sympy as sp
from fpylll import IntegerMatrix,BKZ
from fractions import Fraction
p=pathlib.Path(__file__).resolve().parent/'CTturyn1';d=np.load(p/'frame_parts.npz');P=d['parts'];pd=d['denominators'];n=36
params=[(Fraction(1,3),Fraction(1),Fraction(3)),(Fraction(1,2),Fraction(1),Fraction(5,2)),(Fraction(1,4),Fraction(1),Fraction(13,4))]
for idx,x in enumerate(params):
 den=int(sp.ilcm(*[(q/Fraction(int(t))).denominator for q,t in zip(x,pd)]));Q=sum(int(q/Fraction(int(t))*den)*M for q,t,M in zip(x,pd,P));B0=np.rint(np.linalg.cholesky(Q/den)*1e8).astype(np.int64);B=IntegerMatrix.from_matrix(B0.tolist());U=IntegerMatrix.identity(n);BKZ.reduction(B,BKZ.Param(block_size=28,max_loops=3),U=U,float_type='long double');U=np.array([[int(U[i,j])for j in range(n)]for i in range(n)],np.int64);R=U@Q@U.T;np.savetxt(p/f'candidate{idx}_Q.txt',R,fmt='%d');np.savetxt(p/f'candidate{idx}_R.txt',np.linalg.cholesky(R).T,fmt='%.17g')
 run=subprocess.run([str(pathlib.Path(__file__).resolve().parent.parent/'bin'/'shell_discovery'),'36',str(p/f'candidate{idx}_Q.txt'),str(p/f'candidate{idx}_R.txt'),str(8*den-1),'100000000',str(p/f'candidate{idx}_bad.txt'),'0'],capture_output=True,text=True);r=json.loads(run.stdout);r.update({'scales':list(map(str,x)),'denominator':den,'density_if_minimum8':float(2**18/3**9/np.prod([float(t)for t in x])**6)})
 if (p/f'candidate{idx}_bad.txt').stat().st_size:
  V=np.loadtxt(p/f'candidate{idx}_bad.txt',dtype=np.int64,ndmin=2)[:,1:]@U;norms=np.einsum('ij,jk,ik->i',V,Q,V);r['minimum_upper']=float(norms.min()/den);r['witness']=V[norms.argmin()].tolist();r['witness_profile']=[str(Fraction(int(v@M@v),int(dd)))for v in [V[norms.argmin()]] for M,dd in zip(P,pd)];np.save(p/f'candidate{idx}_bad.npy',V)
 print(json.dumps(r),flush=True);(p/f'candidate{idx}_report.json').write_text(json.dumps(r,indent=2));np.savez(p/f'candidate{idx}.npz',Q=Q,den=den,U=U)
