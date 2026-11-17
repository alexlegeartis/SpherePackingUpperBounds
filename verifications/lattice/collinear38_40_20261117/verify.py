"""Exact independent audit; full predecessor replay is available via --full."""
import gzip,hashlib,json,math,sys,subprocess
from pathlib import Path
from fractions import Fraction as F
import sympy as s
if not __debug__: raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
oldp=p.parent/'rank_two_20261117'
rho=F(1,10**60)
assert sum((F(168)**j/F(math.factorial(j)) for j in range(401)),F(0))>1/rho
Gamma=(1-rho)*(1-rho*rho)
f=oldp/'rank2_candidate42_B115950994_20261117.json'
assert hashlib.sha256(f.read_bytes()).hexdigest()=='3d9ee200d60ab7784a169bf73059a4e7e0688423060f28987121d1828d2162cf'
x=s.symbols('x');raw=json.loads(f.read_text())
Q=s.Poly(sum(s.Rational(str(v))*s.gegenbauer(2*j,20,x)/s.gegenbauer(2*j,20,1)for j,v in enumerate(raw['even_coefficients']))-s.Rational(1,5000),x)
rows=[]
for n,safe,prev,external,lower in [(38,'0.000012713806752','0.00001271981515','0.0000152662254','0.0000002175667308912894917825935149131191'),(39,'0.000008757804412','0.000008763069958','0.0000109486830','0.0000001680281773654323233692192288852917'),(40,'0.000006157949542','0.000006159890437','0.00000782816106','0.0000001479833993134267670269835311956182')]:
 rem=Q;co=[]
 for k in range(Q.degree(),-1,-1):
  b=s.Poly(s.gegenbauer(k,s.Rational(n-2,2),x)/s.gegenbauer(k,s.Rational(n-2,2),1),x)
  a=rem.nth(k)/b.nth(k);assert a>=0;co.append(a);rem-=a*b
 assert rem.is_zero and co[-1]>0
 q=F(Q.eval(1)/co[-1])
 blob=gzip.decompress((oldp/f'certificate{n}.json.gz').read_bytes());old=json.loads(blob)
 assert hashlib.sha256(blob).hexdigest()==json.loads((oldp/'MANIFEST.json').read_text())[f'certificate{n}.json']
 ang=old['cyclic_angular'];a=list(map(F,ang['gegenbauer_coefficients']));B=sum(a);assert B==F(ang['B_exact'])
 assert F(ang['profile_a_exact'])==-96 and F(ang['profile_b_exact'])==120
 c=F(old['C'])/B**2;u0=c*sum(a[::2]);v0=c*F(old.get('origin_ratio','1'))*B
 origin0=F(old.get('origin_C',old['C']));assert origin0==v0*v0/c
 u=c*B*Gamma/2;v=v0+u-u0
 assert 0<u0<=u and 2*u<=c*B*Gamma
 origin=q*v*v/u
 extra=2*(2**n-1)*F(old['kappa'])+F(old['extra_objective'])
 obj=F(old['p'][0])+origin0+extra;newobj=F(old['p'][0])+origin+extra
 assert 0<newobj<obj
 bound=F(old['safe_density'])*newobj/obj
 assert bound<F(safe)<F(prev)<F(external)
 rows.append(dict(dimension=n,upper=safe,exact_upper=str(bound),previous_registered_upper=prev,external_upper=external,external_source='Zubrilina (2026), Table 1, https://arxiv.org/html/2610.35830v1#S4.T1',lower=lower,internal_improvement_percent=float(100*(1-F(safe)/F(prev))),external_improvement_percent=float(100*(1-F(safe)/F(external))),external_gap_closed_percent=float(100*(F(external)-F(safe))/(F(external)-F(lower))),u=str(u),v=str(v),origin=str(origin),prior_sha256=hashlib.sha256(blob).hexdigest()))
 print(n,safe,'EXACT PASS')
(p/'RESULTS.json').write_text(json.dumps(rows,indent=2)+'\n')
if '--full' in sys.argv:
 for n in (38,39,40):subprocess.run([sys.executable,str(oldp/'run_all.py'),'--dimension',str(n)],check=True)
print('COLLINEAR TRANSFER PASS')
