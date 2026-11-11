"""Verify the frozen unrestricted five-dimensional density certificate.
Run without Python optimization flags. No solver or network is required.
"""
if not __debug__: raise RuntimeError('Verification must not run under python -O.')
import json,sys,subprocess,hashlib
from pathlib import Path
from fractions import Fraction
from flint import fmpq as Q,fmpq_mat as Mat,fmpq_poly as Poly,arb,ctx
from gaussian_guards import verify_guards,PI_LO
ctx.prec=256
ROOT=Path(__file__).resolve().parent
s=json.loads((ROOT/'witness.json').read_text());z=json.loads((ROOT/'witness_psd_exact.json').read_text());pre=json.loads((ROOT/'prerequisites.json').read_text())
def fq(v): return Q(str(Fraction(str(v))))
assert hashlib.sha256((ROOT/'source_exact.json').read_bytes()).hexdigest()=='235920216afbd10b38be5efd3672e9d94b5ec5af8fc54aca3a63d77d9997721f'
raw=json.loads((ROOT/'source_exact.json').read_text());assert len(raw)==1;raw=raw[0]
assert (raw['dimension'],raw['radial_degree'],raw['radial_nodes'],raw['angular_degree'])==(5,31,5,12)
assert fq(raw['exclusion_squared'])==1 and fq(raw['truncation_squared'])==Q(9,4)
assert z['source']=='witness.json' and z['matrices']==raw['matrices']
assert s['x']==raw['coefficients'] and s['count_weights']==raw['count_weights'] and s['count_bounds']==raw['count_bounds']
assert pre['dimension']==5 and s['lens_budget']=='step_allowances_005.json'
g=pre['guards']['reciprocal_gaussian_guards']
assert g['eta']==g['zeta']=='0' and g['physical_intervals']==g['fourier_intervals']==[]
assert fq(g['c'])==Q(157,50)
budget=json.loads((ROOT/'step_allowances_005.json').read_text())
assert budget['schema']=='tail_funded_lens_step_allowances_v1'
records=[r for r in budget['dimensions'] if r['dimension']==5];assert len(records)==1
record=records[0]
assert record.get('geometry_mode','centroid')=='centroid' and fq(record['short_radius'])==Q(27,25)
assert fq(record['width'])==Q(1,200) and len(record['bins'])==104
assert [(fq(t['lower_base']),t['cap']) for t in record['exact_cap_geometry']]==[(Q(41,25),83),(Q(17,10),8),(Q(7,4),5),(Q(9,5),3),(Q(46,25),2),(Q(48,25),1)]
assert all(fq(t['hi'])-fq(t['lo'])==Q(1,200) for t in record['bins'])
payload=record['physical_function']
assert payload==dict(dimension=5,basis='L_k^(n/2-1)(2X)/L_k^(n/2-1)(0)',variable='X=pi*r^2',physical_factor='exp(-X)',additional_physical_terms='none',degree=31,exact_rational_coefficients=[str(Fraction(str(c))) for c in s['x']])
from step_allowance_query import self_test
assert self_test(record)
print('Frozen source, exact payload, table and boundary conventions: PASS',flush=True)

assert s['dimension']==5 and s['J']==1 and s['p']==4 and s['K']==12
assert s['kernel_family']=='tail_funded_lens_cyclic' and fq(s['R'])==Q(3,2)
assert list(map(fq,s['count_radii']))==list(map(fq,['1.003012','1.008096','1.013258','1.017445']))
assert list(map(fq,s['count_bounds']))==list(map(Q,[45,46,47,48]))
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
X=Poly([0,1]);lag=[Poly([1]),Q(5,2)-2*X]
for k in range(1,31):lag.append(((2*k+Q(5,2)-2*X)*lag[-1]-(k+Q(3,2))*lag[-2])/Q(k+1))
P=sum((coeff[k]*lag[k]/lag[k][0]for k in range(32)),Poly([]));H=sum(((-1)**k*coeff[k]*lag[k]/lag[k][0]for k in range(32)),Poly([]))
assert pre['guards']['dimension']==5
addition,normalization=verify_guards(pre['guards'],P,H,Q(1),Q(9,4),PI_LO)
B=Mat([[M[0][i+1,j+1]for j in range(5)]for i in range(5)]);b=Mat([[M[0][i+1,0]]for i in range(5)]);v=B.solve(b);sat=(b.transpose()*v)[0,0];alpha=M[0][0,0]-sat;assert alpha>=0 and alpha==Q(pre['origin_slack'])
W=Mat([[v[j,0]for j in range(5)]]+[[Q(int(i==j))for j in range(5)]for i in range(5)])
S=W*B*W.transpose();assert S[0,0]==sat
assert all(S[i,j]==M[0][i,j]for i in range(6)for j in range(6)if (i,j)!=(0,0))
count_cost=sum((l*Q(k)for l,k in zip(lambdas,[45,46,47,48])),Q(0))
objective=(P[0]+sat+count_cost+addition)/normalization
assert objective==Q(pre['saturated_objective'])
density=arb.pi()**2/60*arb(objective)
assert density<arb('0.509996147')
print('Exact Schur factorization and outward objective: PASS',density,flush=True)
subprocess.run([sys.executable,str(ROOT/'audit_algebra_exact.py')],check=True)
subprocess.run([sys.executable,str(ROOT/'construct_step_allowances.py'),'--verify',str(ROOT/'step_allowances_005.json'),'--candidate',str(ROOT/'source_exact.json')],check=True)
subprocess.run([sys.executable,str(ROOT/'verify_count_endpoints.py')],check=True)
for stem,suffix in [('certify_count_pair_arb','count_pair_arb'),('certify_cyclic_bernstein','cyclic_bernstein')]:
 subprocess.run([sys.executable,str(ROOT/f'{stem}.py'),str(ROOT/'witness_psd_exact.json')],check=True)
 check=json.loads((ROOT/f'witness_{suffix}.json').read_text());assert check['verified'],(stem,check)
subprocess.run([sys.executable,str(ROOT/'certify_wedge_far.py'),str(ROOT/'witness_psd_exact.json'),str(ROOT/'step_allowances_005.json')],check=True)
far=json.loads((ROOT/'witness_wedge_far.json').read_text());assert far['verified']
out={'dimension':5,'verified':True,'strict_upper':'0.509996147','density_enclosure':str(density),'scope':'Unrestricted sphere packings; published near-neighbor count bounds are external mathematical inputs. Elementary centroid lens caps and tail budgets are replayed.'}
(ROOT/'verification_result.json').write_text(json.dumps(out,indent=2)+'\n');print('COMPLETE_CERTIFICATE_PASS',json.dumps(out),flush=True)
