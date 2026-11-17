from pathlib import Path
import sympy as S,json,hashlib
p=Path(__file__).parent;f=p/'input/frozen.json';z=json.loads(f.read_text());x=S.symbols('x');c=list(map(S.Rational,z['co']));B=[S.Matrix([[S.Rational(t)for t in row]for row in b])for b in z['B']]
# Explicit quadratic Bernstein controls, independently of producer matrix transforms.
vectors=set()
for j in range(8):
 a=S.Rational(j,8);b=a+S.Rational(1,8)
 vectors.update([(1,a,a*a),(1,(a+b)/2,a*b),(1,b,b*b)])
polys=[S.Poly(S.gegenbauer(2*k,S.Rational(20),x)/S.gegenbauer(2*k,S.Rational(20),1),x)for k in range(len(c))]
count=0
for i,v in enumerate(sorted(vectors)):
 for w in sorted(vectors)[i:]:
  poly=sum((polys[k]*c[k]*(S.Matrix(v).T*B[k]*S.Matrix(w))[0]for k in range(len(c))),S.Poly(0,x))+S.Poly(S.Rational(1,10**7),x)
  assert poly.eval(0)<0 and poly.count_roots(0,S.Rational(1,2))==0
  count+=1
assert len(vectors)==17 and count==153
out={'status':'INDEPENDENT EXACT CAP PASS','polynomials':count,'controls':len(vectors),'cosine_interval':'[0,1/2]','margin':'1/10000000','frozen_sha256':hashlib.sha256(f.read_bytes()).hexdigest()};(p/'cap_independent_result.json').write_text(json.dumps(out,indent=2));print(out)
