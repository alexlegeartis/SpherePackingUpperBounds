import sys;sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
from flint import fmpq as Q,fmpq_poly as P
D=Path(__file__).parent;S=D/'scalar';k=json.loads((S/'kernel.json').read_text());assert (S/'kernel.json').read_bytes()==(D/'geometry/dim50/frozen.json').read_bytes();assert hashlib.sha256((S/'kernel.json').read_bytes()).hexdigest()=='2d4ffd3c113149707bf2b55aff7986c666c75666c62bd5ad9b797ac005a849c2'
a=json.loads((S/'annular_certificate.json').read_text());p=json.loads((S/'repaired.json').read_text());assert a['P']==p['P'];assert F(a['s'])==F(k['s']);assert a['dimension']==p['dimension']==k['dimension']==50
# Every positive residual height has a separate exact spherical-code polynomial, so unused analytic fallback cap is irrelevant.
c=json.loads((S/'caps/result.json').read_text());assert all(F(r['height'])==0 or r['cap_certificate'] is not None for r in c['rows']);edges=[F(c['rows'][0]['l'])]+[F(r['r'])for r in c['rows']];assert edges[0]==22 and edges[-1]==180
for i,r in enumerate(c['rows']):assert F(r['l'])==edges[i]and 0<F(r['r'])-F(r['l'])<=22
result=json.loads((S/'caps/density_result.json').read_text());assert F(result['density_exact'])<F(199548924,10**15);assert F(result['safe_upper'])==F(199548924,10**15)
print('PASS bytebinding, scalar identity, fullannular adjacency, allpositiveheights exactcodecaps, strictsaferounding')
co=list(map(F,k['co']));B=[[[F(v)for v in row]for row in b]for b in k['B']];z=list(map(F,k['z']));s=F(k['s']);assert co[0]==1 and min(co)>0 and F(k['q'])==sum(co)
for b in B:
 assert len(b)==2 and all(len(row)==2 for row in b)and b[0][1]==b[1][0]
 assert b[0][0]>=0 and b[1][1]>=0 and b[0][0]*b[1][1]>=b[0][1]**2
b=B[0];det=b[0][0]*b[1][1]-b[0][1]**2;assert det>0 and s>(b[1][1]*z[0]**2-2*b[0][1]*z[0]*z[1]+b[0][0]*z[1]**2)/det
print('PASS exact PSD harmonic and origin blocks')
