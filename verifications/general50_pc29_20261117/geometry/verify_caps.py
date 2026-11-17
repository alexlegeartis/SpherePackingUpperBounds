import json,sys,hashlib
from pathlib import Path
from fractions import Fraction as F
from flint import fmpq as Q,fmpq_poly as P,arb,arb_poly,ctx
if not __debug__:raise RuntimeError('Assertions required')
sys.set_int_max_str_digits(0);ctx.prec=512;n=int(sys.argv[1]);D=Path(__file__).parent/f'dim{n}';r=json.load(open(D/'frozen.json'));co=list(map(Q,r['co']));B=[[[Q(v)for v in row]for row in b]for b in r['B']];x=P([0,1]);G=[P([1]),x]
for k in range(2,len(co)):G.append(((2*k+n-4)*x*G[-1]-(k-1)*G[-2])/(k+n-3))
M=[[sum((co[k]*B[k][i][j]*G[k]for k in range(len(co))),P([0]))for j in range(2)]for i in range(2)];out=[]
for p in [M[0][0],M[0][0]+M[0][1],M[0][0]+2*M[0][1]+M[1][1]]:
 ap=arb_poly([arb(v)for v in p]);vals=[arb(p(-1)),arb(p(Q(1,2)))];roots=p.derivative().complex_roots();assert sum(m for z,m in roots)==p.degree()-1
 for z,m in roots:
  if not z.imag.contains(0)or z.real < -1 or z.real > arb(1)/2:continue
  vals.append(ap(z.real))
 assert all(v<0 for v in vals),list(map(str,vals));out.append(dict(polynomial=list(map(str,p)),roots_total=sum(m for z,m in roots),tested=len(vals),maximum_upper=str(max(v.upper()for v in vals))))
(D/'caps_result.json').write_text(json.dumps(dict(status='FULL[-1,1/2] CORNER SIGNS PASS',precision=512,frozen_sha=hashlib.sha256((D/'frozen.json').read_bytes()).hexdigest(),corners=out),indent=2));print([v['maximum_upper']for v in out])
