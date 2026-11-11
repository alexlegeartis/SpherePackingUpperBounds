"""Verify the frozen unrestricted six-dimensional density certificate.
Run without Python optimization flags. No solver or network is required.
"""
if not __debug__: raise RuntimeError('Verification must not run under python -O.')
import json,sys,subprocess
from pathlib import Path
from fractions import Fraction
from flint import fmpq as Q,fmpq_mat as Mat,fmpq_poly as Poly,arb,ctx
from gaussian_guards import verify_guards,PI_LO
ctx.prec=256
ROOT=Path(__file__).resolve().parent
s=json.loads((ROOT/'witness.json').read_text());z=json.loads((ROOT/'witness_psd_exact.json').read_text());pre=json.loads((ROOT/'prerequisites.json').read_text())
def fq(v): return Q(str(Fraction(str(v))))
assert s['dimension']==6 and s['J']==1 and s['p']==4 and s['K']==12
assert s['kernel_family']=='chebyshev_radius_count_cyclic' and fq(s['R'])==Q(3,2)
assert list(map(fq,s['count_radii']))==list(map(fq,['1.002005','1.004023','1.006053','1.008096']))
assert list(map(fq,s['count_bounds']))==list(map(Q,[79,80,81,82]))
lambdas=list(map(fq,s['count_weights']));assert len(lambdas)==4 and min(lambdas)>=0
coeff=list(map(fq,s['x']));assert len(coeff)==32
M=[Mat([[Q(v)for v in row]for row in block])for block in z['matrices']];assert len(M)==13
for ell,A in enumerate(M):
 N=6 if ell==0 else 5;assert A.nrows()==A.ncols()==N and A==A.transpose()
 L=Mat(N,N);diag=[]
 for j in range(N):
  L[j,j]=1;d=A[j,j]-sum((L[j,k]**2*diag[k]for k in range(j)),Q(0));assert d>0;diag.append(d)
  for i in range(j+1,N):L[i,j]=(A[i,j]-sum((L[i,k]*L[j,k]*diag[k]for k in range(j)),Q(0)))/d
print('Exact rational positive definiteness: PASS',flush=True)
X=Poly([0,1]);lag=[Poly([1]),3-2*X]
for k in range(1,31):lag.append(((2*k+3-2*X)*lag[-1]-(k+2)*lag[-2])/Q(k+1))
P=sum((coeff[k]*lag[k]/lag[k][0]for k in range(32)),Poly([]));H=sum(((-1)**k*coeff[k]*lag[k]/lag[k][0]for k in range(32)),Poly([]))
assert pre['guards']['dimension']==6
addition,normalization=verify_guards(pre['guards'],P,H,Q(1),Q(9,4),PI_LO)
B=Mat([[M[0][i+1,j+1]for j in range(5)]for i in range(5)]);b=Mat([[M[0][i+1,0]]for i in range(5)]);v=B.solve(b);sat=(b.transpose()*v)[0,0];alpha=M[0][0,0]-sat;assert alpha>=0 and alpha==Q(pre['origin_slack'])
W=Mat([[v[j,0]for j in range(5)]]+[[Q(int(i==j))for j in range(5)]for i in range(5)])
S=W*B*W.transpose();assert S[0,0]==sat
assert all(S[i,j]==M[0][i,j]for i in range(6)for j in range(6)if (i,j)!=(0,0))
count_cost=sum((l*Q(k)for l,k in zip(lambdas,[79,80,81,82])),Q(0))
objective=(P[0]+sat+count_cost+addition)/normalization
assert objective==Q(pre['saturated_objective'])
density=arb.pi()**3/384*arb(objective)
assert density<arb('0.409128533')
print('Exact Schur factorization and outward objective: PASS',density,flush=True)
for stem,suffix in [('certify_count_pair_arb','count_pair_arb'),('certify_far_bernstein','far_bernstein'),('certify_cyclic_bernstein','cyclic_bernstein')]:
 subprocess.run([sys.executable,str(ROOT/f'{stem}.py'),str(ROOT/'witness_psd_exact.json')],check=True)
 check=json.loads((ROOT/f'witness_{suffix}.json').read_text());assert check['verified'],(stem,check)
out={'dimension':6,'verified':True,'strict_upper':'0.409128533','density_enclosure':str(density),'scope':'Unrestricted sphere packings; published near-neighbor count bounds are external mathematical inputs.'}
(ROOT/'verification_result.json').write_text(json.dumps(out,indent=2)+'\n');print('COMPLETE_CERTIFICATE_PASS',json.dumps(out),flush=True)
