import json,math,time,hashlib,sympy as s
from pathlib import Path
D=Path(__file__).parent;z=json.load(open(D/'frozen.json'));co=list(map(s.Rational,z['co']));Bs=[s.Matrix([[s.Rational(v)for v in row]for row in M])for M in z['B']];x=s.symbols('x');u=s.symbols('u');P=[]
for i in range(3):
 row=[]
 for j in range(3):
  pp=s.Poly(sum(co[k]*Bs[k][i,j]*s.gegenbauer(2*k,s.Rational(41,2),x)/s.gegenbauer(2*k,s.Rational(41,2),1)for k in range(len(co))),x);assert all(pp.nth(k)==0 for k in range(1,pp.degree()+1,2));row.append(s.Poly(sum(pp.nth(2*k)*u**k for k in range(len(co))),u))
 P.append(row)
L=s.Matrix([[s.Rational(math.comb(i,j),math.comb(2,j))if j<=i else 0 for j in range(3)]for i in range(3)]);vec=set()
for a in range(8):
 l=s.Rational(a,8);A=s.Matrix([[math.comb(k,j)*l**(k-j)*s.Rational(1,8)**j if k>=j else 0 for k in range(3)]for j in range(3)]);T=L*A
 for i in range(3):vec.add(tuple(T[i,j]for j in range(3)))
vec=sorted(vec);start=time.time();count=0;eta=s.Rational(1,10**7)
for ii,v in enumerate(vec):
 for jj in range(ii,len(vec)):
  w=vec[jj];poly=sum((P[i][j]*v[i]*w[j]for i in range(3)for j in range(3)),s.Poly(0,u))+s.Poly(eta,u)
  assert poly.eval(0)<0 and poly.count_roots(0,s.Rational(1,4))==0,(ii,jj,str(poly.eval(0)))
  count+=1
 print(ii,count,time.time()-start,flush=True)
out=dict(status='EXACT CAP BERNSTEIN-STURM PASS',control_vectors=len(vec),angular_polynomials=count,h_subintervals=8,angular_square_interval=['0','1/4'],negative_margin=str(eta),input_sha=hashlib.sha256((D/'frozen.json').read_bytes()).hexdigest());(D/'cap_exact_result.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
