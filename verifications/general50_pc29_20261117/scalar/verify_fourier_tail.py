from pathlib import Path
import json,math,sys
from flint import fmpq as Q,fmpq_poly as Poly,ctx
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);ctx.prec=1024;D=Path(__file__).parent
r=json.loads((D/'repaired.json').read_text());assert r['dimension']==50 and r['contact_x']=='22'
p=Poly(list(map(Q,r['P'])));q=Poly(list(map(Q,r['Q'])));factorized=Poly([1])
for f in r['factors']:
 if f['type']=='negative_linear':
  a=Q(f['a']);assert a>0;factorized*=Poly([1,1/a])
 elif f['type']=='positive_double':
  a=Q(f['mid']);assert a>0;factorized*=Poly([1,-2/a,1/(a*a)])
 elif f['type']=='complex_pair':
  re,im=Q(f['re']),Q(f['im']);assert im>0;norm=re*re+im*im;factorized*=Poly([1,-2*re/norm,1/norm])
 else:raise AssertionError('Unknown factor')
assert factorized==q and q(0)==1 and p.degree()==q.degree()==71
reconstructed=[]
for j in range(p.degree()+1):
 value=Q(0)
 for k in range(j,p.degree()+1):
  rising=Q(1)
  for ell in range(k-j):rising*=Q(50,2)+j+ell
  value+=p[k]*math.factorial(k)*(-1)**j*rising/(math.factorial(k-j)*math.factorial(j))
 reconstructed.append(value)
assert Poly(reconstructed)==q
print('EXACT FACTORS AND SEPARATE CLOSED-FORM FOURIER CHECK PASS',flush=True)
roots=p.complex_roots();assert sum(m for z,m in roots)==p.degree()
assert all(not z.imag.contains(0) or z.real<180 for z,m in roots)
assert p(Q(180))<0
out=dict(status='SEPARATE EXACT FACTORS, TRANSFORM, FULL SPATIAL TAIL PASS',precision=1024,dimension=50,degree=71,tail_cut_x=180,total_root_multiplicity=sum(m for z,m in roots),possibly_real_root_multiplicity=sum(m for z,m in roots if z.imag.contains(0)),all_possible_real_roots_below_cut=True)
(D/'tail_transform_result.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
