from pathlib import Path
from fractions import Fraction as F
import json,hashlib,math
D=Path(__file__).parent
out=[]
for n in [44,45]:
 p=D/f'input{n}';a=json.load(open(p/f'root_repair{n}/input_candidate.json'));g=json.load(open(p/'frozen.json'));q=F(g['q']);assert F(a['q'])==q and F(a['a'])>0
 C=[[[F(v)for v in row]for row in M]for M in a['C']]
 if n==44:
  co=list(map(F,g['co']));B=[[[F(v)for v in row]for row in M]for M in g['B']];assert C==[[[co[k]*v for v in row]for row in M]for k,M in enumerate(B)];assert list(map(F,a['z']))==list(map(F,g['z'])) and F(a['s'])==F(g['s'])
 else:
  r=json.load(open(p/'restriction45.json'));assert C==[[[F(v)for v in row]for row in M]for M in r['C']];assert r['source_sha']==hashlib.sha256((p/'frozen.json').read_bytes()).hexdigest()
 z=list(map(F,a['z']));S=F(a['s']);B=[[F(a['a'])*v for v in row]for row in C[0]];det=B[0][0]*B[1][1]-B[0][1]**2;assert det>0;schur=S-(B[1][1]*z[0]**2-2*B[0][1]*z[0]*z[1]+B[0][0]*z[1]**2)/det;assert schur>0
 for M in C:assert M[0][1]==M[1][0] and M[0][0]>=0 and M[1][1]>=0 and M[0][0]*M[1][1]>=M[0][1]**2
 ann=json.load(open(p/f'root_repair{n}/annular_final_certificate.json'));rep=json.load(open(p/f'root_repair{n}/repaired.json'));assert ann['P']==rep['P'];assert ann['C']==a['C'] and ann['z']==a['z'] and F(ann['s'])==S
 last=F(22);budget=F(0)
 for row in ann['rows']:
  l,r,h=F(row['l']),F(row['r']),F(row['height']);assert l==last and r>l and r-l<=88 and h>=0;last=r;R=r/22
  if R<2:cap=1
  else:
   t=1-2/R;mu=F(1);caps=[]
   for k in range(501):
    rat=F(2*k+1,n+2*k);den=mu*(rat-t*t)
    if den>0:caps.append((1-t*t)/den)
    mu*=rat
   cap=int(min(caps))
  assert cap==row['cap'];budget+=2**(n+1)*cap*h
 assert last==180 and budget==F(ann['budget'])
 print(n,'exact geometry/scalar links PASS, Schur',schur);out.append(dict(n=n,schur=str(schur),q=str(q)))
(D/'link_result.json').write_text(json.dumps(out,indent=2))
