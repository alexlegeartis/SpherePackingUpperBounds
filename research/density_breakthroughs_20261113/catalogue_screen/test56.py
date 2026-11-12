#!/usr/bin/env python3
"""Search for short-vector witnesses along the exact 56D deformation.

The returned witnesses are checked exactly. Failure to find one is not proof.
"""
from pathlib import Path
from fractions import Fraction
import sys,json,time,subprocess
import numpy as np
from fpylll import IntegerMatrix,BKZ,GSO,LLL

P=Path(__file__).resolve().parent
d=np.load(P/'L56a_flex.npz');G=d['G'];D=d['D'];n=len(G)
t=Fraction(sys.argv[1] if len(sys.argv)>1 else '3');cap=int(sys.argv[2])if len(sys.argv)>2 else 300000000
Q=t.denominator*G+t.numerator*D
assert np.linalg.eigvalsh(Q)[0]>0
start=time.time();B0=np.rint(np.linalg.cholesky(Q)*1e7).astype(np.int64)
B=IntegerMatrix.from_matrix(B0.tolist());U=IntegerMatrix.identity(n)
BKZ.reduction(B,BKZ.Param(block_size=28,max_loops=3),U=U,float_type='long double')
U=np.array([[int(U[i,j])for j in range(n)]for i in range(n)],dtype=np.int64);Qr=U@Q@U.T
q=np.diag(Qr);print('t',str(t),'shortest_reduced_basis_norm',float(q.min()/t.denominator),'seconds',time.time()-start,flush=True)
tag=str(t).replace('/','_').replace('-','m');out=P/f'test56_{tag}';out.mkdir(exist_ok=True)
np.savetxt(out/'Q.txt',Qr,fmt='%d');np.savetxt(out/'R.txt',np.linalg.cholesky(Qr).T,fmt='%.17g');np.savez_compressed(out/'reduction.npz',Q=Q,Qr=Qr,U=U)
if q.min()<28*t.denominator:
    witnesses=U[q<28*t.denominator]
    enum={'skipped':'short basis vector already rejects candidate'}
else:
    run=subprocess.run([str(P.parent/'bin/shell_discovery'),str(n),str(out/'Q.txt'),str(out/'R.txt'),str(28*t.denominator-1),str(cap),str(out/'bad.txt'),'1'],capture_output=True,text=True)
    enum=json.loads(run.stdout)
    witnesses=np.loadtxt(out/'bad.txt',dtype=np.int64,ndmin=2)[:,1:]@U if (out/'bad.txt').stat().st_size else np.empty((0,n),dtype=np.int64)
rows=[]
for v in witnesses:
    a=int(v@G@v);b=int(v@D@v);assert Fraction(a)+t*b<28
    rows.append({'vector':v.tolist(),'seed_norm':a,'slope':b,'candidate_norm':str(Fraction(a)+t*b),'barrier_t':str(Fraction(28-a,b)) if b else None})
report={'t':str(t),'enumeration':enum,'witnesses':rows,'conditional_center_density':float(7**28/(29-float(t)**2)**14),'seconds':time.time()-start,'scope':'discovery; every listed witness exact; no empty-search minimum certification'}
(out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
