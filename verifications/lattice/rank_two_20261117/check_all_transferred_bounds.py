import json,hashlib,sympy as s
from fractions import Fraction as F
from pathlib import Path
p=Path(__file__).parent;x=s.symbols('x');rho=F(1,100);Gamma=(1-rho)*(1-rho*rho)
files={43:'rank2_candidate43_B171069097_20261117.json',42:'rank2_candidate42_B115950994_20261117.json',45:'rank2_candidate45_B318469503_20261117.json',48:'rank2_candidate_B783489561_20261117.json'}
expected={43:'c2481b48633f11417291c75bb6a44c833fa4ea8d1277a3de1fa2e4e79ccf195b',42:'3d9ee200d60ab7784a169bf73059a4e7e0688423060f28987121d1828d2162cf',45:'defe0aa8ab8d91f7bd14a983c341bf749c3695d01af509431c5ae77e0b6f1712',48:'119dc3c5ef12da19e6edb01332e87e4c472da94d6c60c23fa438f74b7da147ae'}
polys={}
for k,file in files.items():
 raw=(p/file).read_bytes();assert hashlib.sha256(raw).hexdigest()==expected[k];co=json.loads(raw)['even_coefficients'];a=s.Rational(k-2,2);Q=s.Poly(sum(s.Rational(str(v))*s.gegenbauer(2*j,a,x)/s.gegenbauer(2*j,a,1)for j,v in enumerate(co))-s.Rational(2,10000),x);polys[k]=Q
rows=[]
for n,k,safe in [(38,42,'0.00001271981515'),(39,42,'0.000008763069958'),(40,42,'0.000006159890437'),(41,42,'0.000004493661175'),(42,42,'0.000003521717941'),(43,43,'0.000002554917283'),(44,45,'0.000001779651377'),(45,45,'0.000001180784470'),(46,45,'0.0000007870493008'),(47,48,'0.0000005258790757'),(48,48,'0.0000003529712022')]:
 Q=polys[k];rest=Q;co=[]
 for j in range(Q.degree()//2,-1,-1):
  basis=s.Poly(s.gegenbauer(2*j,s.Rational(n-2,2),x)/s.gegenbauer(2*j,s.Rational(n-2,2),1),x);a=rest.nth(2*j)/basis.nth(2*j);assert a>=0;co.append(a);rest-=a*basis
 assert rest.is_zero and co[-1]>0;qn=F(Q.eval(1)/co[-1])
 raw=(p/f'certificate{n}.json').read_bytes();old=json.loads(raw);ang=old['cyclic_angular'];a=list(map(F,ang['gegenbauer_coefficients']));B=sum(a);assert B==F(ang['B_exact']);qold=sum(a[::2]);c=F(old['C'])/B**2;R=F(old.get('origin_ratio','1'))*B;u0=c*qold;v0=c*R
 assert F(old.get('origin_C',old['C']))==c*R**2 and F(ang['profile_a_exact'])==-96 and F(ang['profile_b_exact'])==120
 u=c*B*Gamma/2 if n in (38,39,40) else u0;v=v0+u-u0;assert u>=u0>0 and 2*u<=c*B*Gamma
 origin=qn*v*v/u;extra=2*(2**n-1)*F(old['kappa'])+F(old['extra_objective']);obj=F(old['p'][0])+F(old.get('origin_C',old['C']))+extra;newobj=F(old['p'][0])+origin+extra;upper=F(old['safe_density'])*newobj/obj
 assert upper<=F(safe)<F(old['safe_density'])
 rows.append(dict(dimension=n,kernel_dimension=k,verified_upper=safe,previous_verified_upper=old['safe_density'],decrease_percent=float(100*(1-newobj/obj)),prior_certificate_sha256=hashlib.sha256(raw).hexdigest(),kernel_sha256=expected[k],new_origin_exact=str(origin)))
 print(n,safe,'PASS',float(100*(1-newobj/obj)))
(p/'VERIFIED_RESULTS.json').write_text(json.dumps(rows,indent=2))
