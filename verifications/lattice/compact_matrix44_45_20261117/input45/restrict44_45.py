from pathlib import Path
import json,sympy as s,hashlib
D=Path(__file__).parent;assert hashlib.sha256((D/'frozen.json').read_bytes()).hexdigest()=='fe60310a72412cc57e47dc0f0825a0aa7806dd80d3c7c0c045888c6925957dab';z=json.load(open(D/'frozen.json'));co=list(map(s.Rational,z['co']));B=[[[s.Rational(v)for v in row]for row in b]for b in z['B']];x=s.symbols('x');P=[[s.Poly(sum(co[k]*B[k][i][j]*s.gegenbauer(2*k,22,x)/s.gegenbauer(2*k,22,1)for k in range(len(co))),x)for j in range(2)]for i in range(2)]
for n in [44,45]:
 C=[[[s.Rational(0)for j in range(2)]for i in range(2)]for k in range(11)]
 for i in range(2):
  for j in range(2):
   rem=P[i][j]
   for k in range(20,-1,-1):
    g=s.Poly(s.gegenbauer(k,s.Rational(n,2)-1,x)/s.gegenbauer(k,s.Rational(n,2)-1,1),x);v=rem.nth(k)/g.nth(k);rem-=v*g
    if k%2:assert v==0
    else:C[k//2][i][j]=v
   assert rem.is_zero
 for M in C:assert M[0][1]==M[1][0] and M[0][0]>=0 and M[1][1]>=0 and M[0][0]*M[1][1]>=M[0][1]**2
 assert C[0][0][0]*C[0][1][1]>C[0][0][1]**2
 out=dict(n=n,q=z['q'],C=[[[str(v)for v in row]for row in M]for M in C],source_sha=hashlib.sha256((D/'frozen.json').read_bytes()).hexdigest(),status='EXACT PSD RESTRICTION PASS; density requires separate scalar proof');(D/f'restriction{n}.json').write_text(json.dumps(out,indent=2));print(n,[[float(v)for v in row]for row in C[0]],flush=True)
