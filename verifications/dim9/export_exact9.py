"""Exact reconstruction of candidate radial inputs; NOT proof by itself."""
import json,sys
from fractions import Fraction as F
from pathlib import Path
from flint import fmpq as Q,fmpq_poly as P
sys.set_int_max_str_digits(0);D=Path(__file__).parent
r=json.loads((D/'refit9.json').read_text());x=P([0,1]);a=Q(9,2);G=[P([1]),a-2*x]
for k in range(1,len(r['radial'])-1):G.append(((2*k+a-2*x)*G[-1]-(k+a-1)*G[-2])/(k+1))
G=[g/g(0) for g in G]
def q(v):return Q(str(F(float(v))))
p=sum((q(c)*g for c,g in zip(r['radial'],G)),P());h=sum(((-1)**k*q(c)*g for k,(c,g) in enumerate(zip(r['radial'],G))),P());N=h(0);p/=N;h/=N;C=q(r['C'])/N;lam=q(r['lambda_count'])/N;eps=Q(1,10**6)
p=(p+eps*(1-4*x/9))/(1-eps);h=(h-eps*(1-4*x/9))/(1-eps);C=C/(1-eps)*Q(1000001,1000000);lam/=1-eps
out=dict(dimension=9,physical=list(map(str,p.coeffs())),fourier=list(map(str,h.coeffs())),C_exact=str(C),lambda_count_exact=str(lam),threshold_v_exact=str(F('1.003012')**2),epsilon_exact=str(eps),s_exact='9/8',certified=False)
(D/'certificate.json').write_text(json.dumps(out,indent=2))
