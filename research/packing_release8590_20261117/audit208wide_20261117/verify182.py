if not __debug__: raise RuntimeError('Run without -O.')
"""Independent Smith/coarea audit of frozen86/87/88 packings; parent replay.
Published MW92 minimum16 inherited. No worldwide priority claim.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
from math import factorial,prod
import numpy as np,sympy as sp,json,hashlib
from sympy.matrices.normalforms import smith_normal_form
from flint import fmpz_mat,fmpq_mat,fq_default_ctx,fq_default_poly_ctx,fmpz_mod_poly_ctx
P=Path(__file__).resolve().parent.parent/'packing_resume_20261117'
src=(P/'verify.py').read_text();exec(src[src.index('def parent():'):src.index('def rogers_comparison_intervals():')]);G,pr=parent();J=sp.Matrix(G.tolist());Ji=J.inv();base=json.load(open(P/'published_simplex_results.json'))
def atan(q,N):return sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(N)),F(0))
pilo=16*atan(5,80)-4*atan(239,81);pihi=16*atan(5,81)-4*atan(239,80)
inputs=[('candidate182_n86_final_frozen.json', '0137adc29fc5cae1bed9ebe6b05c2e2bb0bbada3e0dcdad9c578425d31cebe1b', '2.84477266080525e-19'), ('candidate182_n87_final_frozen.json', '733f3a7a89f9ba459918b5bcf98c2dc2e1171815baf2a459fd054917460f5178', '2.58037767092447e-19'), ('candidate182_n88_frozen.json', 'e8032155aef7c6145651c847e8ed4ede068b0c4aec4b68f525d65f160ecf9db0', '2.42250297455694e-19'), ('candidate182_n90_frozen.json', '9fe591f351043f78976d92d8b8f8f2e071a354a31835e2b98799b47827d029ba', '1.61584267273781e-19')]
reports=[]
for fn,sha,safe in inputs:
 raw=(P/fn).read_bytes();assert hashlib.sha256(raw).hexdigest()==sha;a=json.loads(raw);B=sp.Matrix(a['normals']);S=sp.Matrix(a['labels']);k=B.rows;n=92-k;m=S.rows
 D=smith_normal_form(B,domain=sp.ZZ);assert all(abs(D[i,i])==1 for i in range(k));assert len(set(map(tuple,S.tolist())))==m
 K=B*Ji*B.T;assert all(K[i,i]==sp.Rational(182,47) for i in range(k));Q=K.inv();pairs=[]
 for i in range(m):
  for j in range(i):
   z=S[i,:]-S[j,:];pairs.append((z*Q*z.T)[0])
 beta=max(pairs);assert beta<16 and str(beta)==a['beta'];detperiod=47**14*K.det();d2=m*m*((16-beta)/4)**n/detperiod;assert str(d2)==a['density_squared'];d=F(str(d2));basefile='published_adaptive_results.json' if n==85 else ('published_mw92_results.json' if n==90 else 'published_simplex_results.json');base=json.load(open(P/basefile));old=F(next(x['density_squared']for x in base['rows']if x['dimension']==n));ratio=d/old
 if n%2==0:vlo=pilo**n/factorial(n//2)**2
 else:vlo=F(2**(n+1))*pilo**(n-1)/prod(range(1,n+1,2))**2
 lo=d*vlo;assert F(safe)**2<lo
 with localcontext()as ctx:
  ctx.prec=70;c=(Decimal(d.numerator)/Decimal(d.denominator)).sqrt();gain=100*((Decimal(ratio.numerator)/Decimal(ratio.denominator)).sqrt()-1)
  reports.append({'dimension':n,'baseline_file':basefile,'file':fn,'sha256':sha,'period_determinant':str(detperiod),'beta':str(beta),'center_density_approx':str(c),'rigorous_occupied_density_lower':safe,'gain_percent_vs_published_approx':str(gain),'all_pair_losses':list(map(str,pairs)),'status':'PASS'})
print(json.dumps({'parent_replay':pr,'rows':reports},indent=2))
